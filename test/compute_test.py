import unittest
import io
import contextlib

from numpy.testing import assert_array_almost_equal
import numpy as np


class CovarianceTest(unittest.TestCase):

    def setUp(self):
        from hvarma.read_input import Data
        self.Z_fn = 'data/B001_Z.sac'
        self.N_fn = 'data/B001_N.sac'
        self.E_fn = 'data/B001_E.sac'
        self.data = Data(Z_fname=self.Z_fn, 
                         N_fname=self.N_fn, 
                         E_fname=self.E_fn)

    def test_crosscovariance1(self):
        from hvarma.read_input import Window
        from hvarma.compute import compute_crosscovariance
        self.set1 = ([ 1.73828596e+08+6.60849194e+07j,  6.32054345e+07+2.07464694e+08j,
                     -2.54828589e+08+3.21324056e+08j, -1.31609498e+08+3.89973541e+07j,
                     -4.28178907e+07+5.46974371e+07j,  9.74915622e+07+4.59257886e+08j,
                      3.34012717e+08+5.02854297e+08j,  3.22029944e+08+3.03034173e+08j],
                    [ 1.73828596e+08+6.60849194e+07j, -4.40177727e+08-8.73210035e+07j,
                     -4.15957206e+08-2.84325980e+08j, -2.06746743e+08-3.46445622e+08j,
                     -1.62372715e+08-7.81128143e+07j, -9.61959690e+06+1.13506659e+08j,
                     -1.57010915e+08+8.09483120e+07j, -3.35587188e+08+1.39089075e+08j])

        self.set1 = (np.array(self.set1[0]), np.array(self.set1[1]))

        win = Window(self.data, 10, 20)
        cov_v_zx, cov_zx_v = compute_crosscovariance(win.dataE, 
                                                    win.dataN, 
                                                    win.dataZ, 
                                                    20, 8)

        assert_array_almost_equal(self.set1[0]/1e8, cov_v_zx/1e8)
        assert_array_almost_equal(self.set1[1]/1e8, cov_zx_v/1e8)

    def test_crosscovariance2(self):
        from hvarma.read_input import Window
        from hvarma.compute import compute_crosscovariance
        self.set2 = ([ 7.68839918e+08+3.20875624e+08j,  9.39971753e+07+3.08121405e+08j,
                     -1.38147329e+08+2.80586218e+08j,  5.16377001e+08+3.15477898e+08j,
                     -4.27215147e+08-4.56944843e+07j,  9.03287709e+07+1.35888032e+08j,
                      5.81748416e+07+1.13503438e+08j, -4.24402531e+08-2.24817347e+08j,
                      1.39887806e+08-2.85326011e+08j, -2.47337333e+08-5.30065816e+08j],
                    [ 7.68839918e+08+3.20875624e+08j, -4.51899785e+08+4.77027879e+07j,
                      2.35583170e+08+3.90563747e+07j,  1.11886283e+08-1.57104037e+08j,
                     -4.26073272e+08-1.69493448e+08j,  2.79147134e+08-1.60749381e+07j,
                     -4.12113849e+08-8.02345606e+07j, -3.42816027e+08+1.12512582e+07j,
                      6.89784484e+07-1.10769146e+08j, -2.94089109e+08-1.79898821e+08j])
        
        self.set2 = (np.array(self.set2[0]), np.array(self.set2[1]))

        win = Window(self.data, 10, 50)
        cov_v_zx, cov_zx_v = compute_crosscovariance(win.dataE, 
                                                    win.dataN, 
                                                    win.dataZ, 
                                                    50, 10)

        assert_array_almost_equal(self.set2[0]/1e8, cov_v_zx/1e8)
        assert_array_almost_equal(self.set2[1]/1e8, cov_zx_v/1e8)

    def test_autocovariance1(self):
        from hvarma.read_input import Window
        from hvarma.compute import compute_autocovariance
        self.set3 = ([ 2.87984429e+09+0.00000000e+00j,  1.85655767e+09+4.10987990e+08j,
                      7.06468559e+08+8.27135659e+08j,  2.42983625e+08+1.07233248e+09j,
                     -1.97140877e+08+9.58876013e+08j, -3.53097863e+08+8.45502744e+08j,
                     -5.48681012e+08+7.22833510e+08j, -1.08682662e+09+5.24977287e+08j],
                    [ 1.03108124e+09,  3.35426972e+08, -1.28336100e+08,  3.81134378e+07,
                     -2.74653371e+07,  1.17148899e+07,  1.08677886e+08, -6.59939340e+07])

        self.set3 = (np.array(self.set3[0]), np.array(self.set3[1]))

        win = Window(self.data, 10, 20)
        cov_x, cov_v = compute_autocovariance(win.dataE, 
                                              win.dataN, 
                                              win.dataZ, 
                                              20, 8)

        assert_array_almost_equal(self.set3[0]/1e8, cov_x/1e8)
        assert_array_almost_equal(self.set3[1]/1e8, cov_v/1e8)

    def test_autocovariance2(self):
        from hvarma.read_input import Window
        from hvarma.compute import compute_autocovariance
        self.set4 = ([ 7.28934951e+09+0.00000000e+00j,  3.26854860e+09+5.90404001e+08j,
                     -3.97073231e+07+1.76561521e+09j, -6.26044627e+07+2.40595417e+09j,
                     -1.44673979e+09+1.78538502e+09j, -1.39247079e+09+7.63443105e+08j,
                     -1.52682997e+09+3.80995884e+08j, -2.37036468e+09+3.75111247e+07j,
                     -1.15803352e+09-4.30143567e+08j,  2.91286747e+08-2.86857765e+08j],
                    [ 2.50563614e+09,  1.73780562e+09,  6.08271785e+08, -5.04970417e+07,
                     -3.43850738e+08, -3.80609258e+08, -7.45510945e+08, -1.22189911e+09,
                     -1.17966176e+09, -6.87825794e+08])
        
        self.set4 = (np.array(self.set4[0]), np.array(self.set4[1]))

        win = Window(self.data, 100, 50)
        cov_x, cov_v = compute_autocovariance(win.dataE, 
                                               win.dataN, 
                                               win.dataZ, 
                                               50, 10)

        assert_array_almost_equal(self.set4[0]/1e8, cov_x/1e8)
        assert_array_almost_equal(self.set4[1]/1e8, cov_v/1e8)


class ModelEquationsTest(unittest.TestCase):

    def setUp(self):
        from hvarma.read_input import Data
        self.Z_fn = 'data/B001_Z.sac'
        self.N_fn = 'data/B001_N.sac'
        self.E_fn = 'data/B001_E.sac'
        self.data = Data(Z_fname=self.Z_fn, 
                         N_fname=self.N_fn, 
                         E_fname=self.E_fn)

    def test_equation_solutions(self):
        from hvarma.read_input import Window
        from hvarma.compute import compute_equations
        self.set1_matrix = np.array([  5.13972248e+19,  3.10707417e+19,  2.11507709e+19, -1.98936107e+18,
                             -1.08894509e+19, -2.35618707e+18, -1.97935396e+18, -4.45278216e+18,
                             -4.21913598e+18, -3.38644522e+18, -9.30082847e+18,  3.10707417e+19,
                              5.13176670e+19,  3.10449962e+19, -5.50313035e+18, -1.95663440e+18,
                             -1.06087876e+19, -2.05685977e+18,  5.00312617e+17, -4.52090459e+18])
        self.set1_indep = np.array([-3.10015175e+19, -2.10644987e+19, -1.53146448e+19,  1.09304118e+19,
                          2.45233710e+18,  2.05290316e+18,  2.07021153e+18, -3.03187250e+18,
                         -2.91581009e+18, -5.90339381e+18, -3.84857288e+18])
        win = Window(self.data, 10, 20)
        mat, indep = compute_equations(win.dataE, win.dataN, win.dataZ, 0.5, 0.5, 100, 3, 50)

        assert_array_almost_equal(self.set1_matrix/1e20, mat.ravel()[:20]/1e20, decimal=8)
        assert_array_almost_equal(self.set1_indep/1e20, indep/1e20, decimal=8)


class ProcessingWindowTest(unittest.TestCase):

    def setUp(self):
        from hvarma.read_input import ArmaParam, Data, Window
        from hvarma.processing import ProcessingWindow
        self.Z_fn = 'data/B001_Z.sac'
        self.N_fn = 'data/B001_N.sac'
        self.E_fn = 'data/B001_E.sac'
        self.filename = 'test/resources/args1.txt'

        self.data = Data(Z_fname=self.Z_fn, 
                         N_fname=self.N_fn, 
                         E_fname=self.E_fn)

        self.param = ArmaParam(filename=self.filename)
        self.window = Window(self.data, 10, self.param.wsize)
        self.pwindow = ProcessingWindow(self.window, self.param)

        self.pwindow.solve_arma() # This should be improved!
        self.pwindow.get_coherence()

    def test_transfer_function(self):
        self.set1 = [0.53337717, 0.49438721, 0.46376495, 0.43970526, 0.42082965, 0.40610029, 0.39473572, 
                    0.38614376, 0.37987166, 0.37556979, 0.37296543, 0.37184367, 0.37203336, 0.37339662, 
                    0.37582087, 0.37921262, 0.3834925, 0.38859114, 0.39444562, 0.40099626, 0.40818365, 
                    0.41594569, 0.4242147, 0.43291446, 0.44195722, 0.45124079, 0.46064582, 0.47003341, 
                    0.47924364, 0.48809509, 0.49638619, 0.50389864, 0.51040347, 0.51566995, 0.51947709, 
                    0.52162723, 0.52196034, 0.52036776, 0.51680332, 0.51129056, 0.50392522, 0.49487308, 
                    0.48436438, 0.47268704, 0.46018151, 0.44724051, 0.43431756, 0.42194862, 0.41079278, 
                    0.40169993]

        transf = self.pwindow.transfer_fun()
        assert_array_almost_equal(transf[:50], self.set1, decimal=5)

    def test_coherence(self):
        self.set2 = [0.31816268, 0.31432628, 0.31126705, 0.3089768, 0.30743961, 0.30663453, 0.30653772, 
                    0.30712409, 0.30836848, 0.3102465, 0.31273501, 0.31581244, 0.31945896, 0.32365648, 0.32838871, 
                    0.33364107, 0.33940065, 0.3456562, 0.35239806, 0.35961814, 0.36730992, 0.37546837, 0.38409, 
                    0.39317273, 0.40271586, 0.41271994, 0.42318655, 0.43411807, 0.44551731, 0.45738693, 0.46972881, 
                    0.48254305, 0.49582677, 0.50957248, 0.52376605, 0.53838407, 0.55339063, 0.56873343, 0.58433904, 
                    0.60010764, 0.6159072, 0.63156763, 0.64687579, 0.66157247, 0.67535328, 0.68787596, 0.69877674, 
                    0.70769826, 0.71432972, 0.71845604]
        coherence = self.pwindow.get_coherence()

        assert_array_almost_equal(coherence[:50], self.set2, decimal=5)

    def test_model_solution(self):

        self.set3_a = [  1.,          -2.72574511,   5.89449761, -10.57235683,  17.90783366,
                       -27.24632772,  39.37763845, -52.8273491,   67.96990399, -82.42774742]
        self.set3_b = [ 0.03996188 +0.17852622j, -0.200071   -0.51689121j,
                      0.56596132 +1.19388204j, -1.12338149 -2.1857838j,
                      1.83750457 +3.67767491j, -2.7198686  -5.65928363j,
                      3.7531891  +8.21177376j, -4.84824045-11.01768525j,
                      5.77761576+14.07635973j, -6.43255569-16.99572242j]

        assert_array_almost_equal(self.pwindow.a[:10], self.set3_a, decimal=6)
        assert_array_almost_equal(self.pwindow.b[:10], self.set3_b, decimal=6)


class AverageDataTest(unittest.TestCase):

    def setUp(self):
        from hvarma.write_output import progress_bar
        from hvarma.read_input import ArmaParam, Data, Window
        from hvarma.processing import ProcessingWindow, AverageData
        
        def get_data_windows(data, size, overlap):
            """ Generator of data slices from data of a given size,
                overlapping one another """
            if data.size < size:
                raise Exception('Window exceeds available data')

            for start in range(0, data.size-size+1, size-overlap):
                yield Window(data, start, size)

        self.Z_fn = 'data/B001_Z.sac'
        self.N_fn = 'data/B001_N.sac'
        self.E_fn = 'data/B001_E.sac'
        self.filename = 'test/resources/args1.txt'
        self.data = Data(Z_fname=self.Z_fn, 
                         N_fname=self.N_fn, 
                         E_fname=self.E_fn)

        param = ArmaParam(filename=self.filename)

        processed_windows = []
        maxwin = 100
        param.maxwin = maxwin
        progress = progress_bar(self.data.size, param.wsize, param.overlap, maxwin)
        f = io.StringIO()
        with contextlib.redirect_stdout(f): # catch stdout
            for idx, data_window in enumerate(get_data_windows(self.data, param.wsize, param.overlap)):
                next(progress)
                model = ProcessingWindow(data_window, param)
                model.solve_arma()
                model.get_coherence()

                processed_windows.append(model)
                if idx+1 == maxwin:  # Limited windows version
                    break

        self.compute_averages = AverageData(processed_windows, param)

    def test_spectrum1(self):
        self.set1 = [0.81721816, 0.8062195,  0.78768426, 0.77677296, 0.77333868, 0.76196565,
                     0.75378387, 0.72431371, 0.70199212, 0.67628692,]

        spectrum = self.compute_averages.get_spectrum_percentile(50)

        assert_array_almost_equal(self.set1, spectrum[:10])


    def test_spectrum2(self):
        self.set2 = [0.52060565, 0.51945786, 0.50757951, 0.5322541,  0.51461029, 0.50894714,
                     0.49031708, 0.47167243, 0.47070937, 0.4549743, ]

        spectrum = self.compute_averages.get_spectrum_percentile(25)

        assert_array_almost_equal(self.set2, spectrum[:10])

    def test_coherence1(self):
        self.set3 = [0.55636983, 0.54936351, 0.54460871, 0.54169427, 0.53825217, 0.53493943,
                     0.53278503, 0.52844035, 0.52371746, 0.51989127]

        coherence = self.compute_averages.get_coherence_percentile(50)

        assert_array_almost_equal(self.set3, coherence[:10])

    def test_coherence2(self):
        self.set4 = [0.47312102, 0.46459376, 0.45645867, 0.44874357, 0.44287254, 0.43955681,
                     0.43276847, 0.42613528, 0.42005411, 0.41519114]

        coherence = self.compute_averages.get_coherence_percentile(25)

        assert_array_almost_equal(self.set4, coherence[:10])

    def test_frequency(self):
        pos_freq, pos_err, neg_freq, neg_err = self.compute_averages.get_frequency(20)

        self.assertAlmostEqual(pos_freq, 3.34310850439)
        self.assertAlmostEqual(pos_err, 0.89931573802)
        self.assertAlmostEqual(neg_freq, -3.34310850439)
        self.assertAlmostEqual(neg_err, 16.539589442)


if __name__ == '__main__':
    unittest.main(verbosity=2)
