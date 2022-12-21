import numpy as np
import ctypes
from numpy.ctypeslib import ndpointer
import os.path
from scipy.linalg import toeplitz


def compute_crosscovariance(dataE, dataN, dataZ, size, maxtau):
    """ Compute crosscovariance according to maximum likelihood estimation. """
    x = dataN + 1j * dataE
    v = dataZ

    # Positive time lags
    cross_cov_zx_v = np.zeros(maxtau, dtype=complex)
    for tau in range(maxtau):
        cross_cov_zx_v[tau] = np.sum(x[tau:] * v[:size - tau]) / size

    # Negative time lags
    cross_cov_v_zx = np.zeros(maxtau, dtype=complex)
    for tau in range(maxtau):
        cross_cov_v_zx[tau] = np.sum(v[tau:] * x[:size - tau]) / size

    return cross_cov_v_zx, cross_cov_zx_v


def compute_autocovariance(dataE, dataN, dataZ, size, maxtau):
    """ Compute autocovariance according to maximum likelihood estimation. """
    x = dataN + 1j * dataE
    v = dataZ

    auto_cov_x = np.zeros(maxtau, dtype=complex)
    for tau in range(maxtau):
        auto_cov_x[tau] = np.sum(x[tau:] * np.conj(x[:size - tau])) / size

    auto_cov_v = np.zeros(maxtau)
    for tau in range(maxtau):
        auto_cov_v[tau] = np.sum(v[tau:] * v[:size - tau]) / size

    return auto_cov_x, auto_cov_v


def compute_equations(dataE, dataN, dataZ, mu, nu, wsize, p, maxtau):
    """ Wrapper of C function to compute equations.
        Uses compiled library "gradient.so".    """
    libname = "./ext_c/gradient.so"
    libpath = os.path.dirname(os.path.abspath(__file__)) + os.path.sep + libname
    lib = ctypes.cdll.LoadLibrary(libpath)
    fun = lib.compute_equations
    fun.restype = None
    fun.argtypes = [ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                    ctypes.c_double, ctypes.c_double,
                    ctypes.c_size_t,
                    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                    ctypes.c_size_t, ctypes.c_int, ctypes.c_int]

    size = 3 * p + 2
    mat = np.zeros((size, size))
    indep = np.zeros(size)

    fun(dataN, dataE, dataZ, mu, nu, size, mat, indep, wsize, p, maxtau)

    return mat, indep


def transfer_function(f0, f1, npun, t, a, b, p):
    """ Compute H/V in frequency range [neg_freq,pos_freq] with freq_points points,
        for a model with coefficients a, b, (length p). t is sampling interval """

    # Check consistency with Nyquist freq.
    if f1 < f0 or 1/(2*t) < max(abs(f0), abs(f1)):
        raise AttributeError('Wrong frequencies')

    freq = np.linspace(f0, f1, npun)
    z = np.exp(-1j * 2 * np.pi * freq * t)
    h = np.zeros(npun, dtype=complex)
    v = np.zeros(npun, dtype=complex)
    for i in range(npun):
        zk = np.power(z[i], np.arange(0, p))
        h[i] = np.dot(zk, b)
        v[i] = np.dot(zk, a)

    return np.abs(h/v)


def compute_coherence(auto_cov_x, auto_cov_v, cross_cov_v_zx, cross_cov_zx_v, nfir, f0, f1, npun, t):
    """ Compute coherence in the [neg_freq, pos_freq] interval.
        nfir is the number of correlation that are considered.
        Read reference for formula details.                     """
    if f1 < f0 or 1 / (2 * t) < max(abs(f0), abs(f1)):
        raise AttributeError('Wrong frequencies')

    # Initialize correlation matrices
    zc1 = toeplitz(auto_cov_x[:nfir])
    zc2 = toeplitz(auto_cov_v[:nfir])
    zsum = zc1 + zc2
    zc12 = toeplitz(cross_cov_zx_v[:nfir], cross_cov_v_zx[:nfir])

    # Perform inversions
    zisum = np.linalg.inv(zsum)
    zic1 = np.linalg.inv(zc1)
    zic2 = np.linalg.inv(zc2)

    # Perform matrix multiplications
    z2isum = np.matmul(zisum, zisum)
    zwork = np.matmul(zisum, zc12)
    znum12 = np.matmul(zwork, zisum)
    z2ic1 = np.matmul(zic1, zic1)
    z2ic2 = np.matmul(zic2, zic2)

    def quadz(a, b):
        """ Quadratic form multiplication of complex vector b with matrix a. """
        return np.conj(b) @ a @ b

    # Compute coherence at each point
    freq = np.linspace(f0, f1, npun)
    z = np.exp(1j * 2 * np.pi * freq * t)
    coh = np.zeros(npun)
    for i in range(npun):
        zste = np.power(z[i], np.arange(0, nfir))  # steering vector
        zcuad1 = quadz(zic1, zste)
        zcuad2 = quadz(zic2, zste)
        zcuad12 = quadz(znum12, zste)
        zcuad11 = quadz(z2ic1, zste)
        zcuad22 = quadz(z2ic2, zste)
        zden = quadz(z2isum, zste)
        coh[i] = np.sqrt(np.abs(zcuad12 / zden) ** 2 /
                         (np.real(zcuad1 / zcuad11) * np.real(zcuad2 / zcuad22)))

    return coh
