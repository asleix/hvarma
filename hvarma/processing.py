"""
Copyright (c) 2022, Spanish National Research Council (CSIC)

Class definitions for processing and results objects.
"""

from functools import lru_cache
from dataclasses import dataclass
from typing import Mapping
import numpy as np
from .compute import compute_crosscovariance, compute_autocovariance,\
                       compute_equations, transfer_function, compute_coherence
from .read_input import ArmaParam, Data


class HVarma:
    """ Class that handles processing of a single time window """
    def __init__(self, window, param):
        """ Initializes window with defined parameters. """
        if not isinstance(window, Data):
            raise AttributeError('Bad data window initialization')
        if not isinstance(param, ArmaParam):
            raise AttributeError('Bad parameter initialization')
        self.data = window.copy()
        self.param = param

        # Center data
        self.muE, self.muN, self.muZ = None, None, None
        self.center()

        # ARMA coefficients
        self.a = None
        self.b = None

        self.coherence = None

    def center(self):
        """ Center series to 0 and assume they are stationary """
        self.muE = np.mean(self.data.dataE)
        self.data.dataE -= self.muE
        self.muN = np.mean(self.data.dataN)
        self.data.dataN -= self.muN
        self.muZ = np.mean(self.data.dataZ)
        self.data.dataZ -= self.muZ

    def solve_arma(self):
        """ Find optimal coefficients for ARMA model minimizing prediction errors. """

        # Find weights
        mu, nu = 0, 0
        if self.param.nu != 'sigma':
            nu = float(self.param.nu)
        if self.param.mu != 'sigma':
            mu = float(self.param.mu)

        # Find system of equations satisfying optimality conditions
        mat, indep = compute_equations(self.data.dataE, self.data.dataN, self.data.dataZ, mu, nu,
                                       self.param.window_size, self.param.model_order, self.param.maxtau)

        solution = np.linalg.solve(mat, indep)

        solution = np.concatenate(([1], solution))

        # Store optimal coefficients
        pp = self.param.model_order+1
        self.a = solution[:pp]
        self.b = solution[pp:2*pp] + 1j*solution[2*pp:]

    def get_correlations(self):
        """ Compute auto and cross correlations of data. """
        auto_cov_x, auto_cov_v = compute_autocovariance(self.data.dataE, self.data.dataN, self.data.dataZ,
                                                        self.param.window_size, self.param.maxtau)
        cross_cov_v_zx, cross_cov_zx_v = compute_crosscovariance(self.data.dataE, self.data.dataN, self.data.dataZ,
                                                                 self.data.size, self.param.maxtau)
        return auto_cov_x, auto_cov_v, cross_cov_v_zx, cross_cov_zx_v

    def transfer_fun(self):
        """ Obtain H/V amplitude from coefficients in the corresponding frequency interval. """
        return transfer_function(self.param.neg_freq, self.param.pos_freq, self.param.freq_points, 1. / self.data.sampling_rate,
                                 self.a, self.b, self.param.model_order + 1)

    def get_coherence(self):
        """ Call corresponding functions to compute coherence. """
        if self.coherence is None:
            auto_cov_x, auto_cov_v, cross_cov_v_zx, cross_cov_zx_v = self.get_correlations()
            self.coherence = compute_coherence(auto_cov_x, auto_cov_v, cross_cov_v_zx, cross_cov_zx_v,
                                               self.param.nfir, self.param.neg_freq, self.param.pos_freq,
                                               self.param.freq_points, 1. / self.data.sampling_rate)
        return self.coherence

    def get_AIC(self):
        """ Compute AIC = n * log (ssr/n) + 2*k. """
        x = self.data.dataN + 1j * self.data.dataE
        v = self.data.dataZ
        p, size = self.param.model_order+1, self.param.window_size
        a, b = self.a, self.b
        ssr = 0
        for i in range(p, size):
            ssr += np.abs(np.dot(x[i - p:i], a) - np.dot(v[i - p:i], b)) ** 2

        AIC = 2 * 3 * p + size * np.log(ssr / size)
        return AIC


class AverageData:
    """ Helper class to handle calculations over all windows """

    def __init__(self, window_list, param):
        assert len(window_list) > 0, "window_list should not be empty"
        self.param = param
        self.num_windows = len(window_list)
        self.station = window_list[0].data.station
        res, coh, aic = [], [], []
        for model in window_list:
            res.append(model.transfer_fun())
            coh.append(model.get_coherence())
            aic.append(model.get_AIC())

        self.spectra = np.vstack(res)
        self.coherence = np.vstack(coh)
        self.AIC = np.array(aic)

    @lru_cache(maxsize=10)
    def get_frequency(self, conf):
        """ Get resonance frequency, corresponding to the maximum peak """
        f0, f1 = self.param.neg_freq, self.param.pos_freq
        freq = np.linspace(f0, f1, self.param.freq_points)

        # Positive peak
        pos_freq, pos_err = 0, 0
        if f1 > 0:
            pos, f = self.spectra[:, freq > 0], freq[freq > 0]
            upp_pos = f[np.argmax(np.percentile(pos, 100-conf/2, axis=0))]
            low_pos = f[np.argmax(np.percentile(pos, conf/2, axis=0))]
            pos_freq = f[np.argmax(np.percentile(pos, 50, axis=0))]
            pos_err = max(abs(upp_pos - pos_freq), abs(pos_freq - low_pos))

        # Negative peak
        neg_freq, neg_err = 0, 0
        if f0 < 0:
            neg, f = self.spectra[:, freq < 0], freq[freq < 0]
            upp_neg = f[np.argmax(np.percentile(neg, 100-conf/2, axis=0))]
            low_neg = f[np.argmax(np.percentile(neg, conf/2, axis=0))]
            neg_freq = f[np.argmax(np.percentile(neg, 50, axis=0))]
            neg_err = max(abs(upp_neg - neg_freq), abs(neg_freq - low_neg))

        return pos_freq, pos_err, neg_freq, neg_err

    def get_spectrum_percentile(self, perc):
        """ Get data corresponding to a given percentile"""
        return np.percentile(self.spectra, perc, axis=0)

    def get_coherence_percentile(self, perc):
        """ Get data corresponding to a given percentile"""
        return np.percentile(self.coherence, perc, axis=0)

    def get_AIC(self):
        """ Return AIC for each window """
        return self.AIC


@dataclass
class OrderSearchResults:
    order_results: Mapping[int, AverageData]
    tol: float
    method: str
    final_order: int
    station: str
    success: bool
