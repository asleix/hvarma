# file to test the module api calls from the end user
import unittest
from numpy.testing import assert_array_almost_equal


class RunModelTest(unittest.TestCase):

    def setUp(self):
        from hvarma import Data, ArmaParam, run_model
        self.param = ArmaParam('test/resources/args1.txt')
        self.data = Data('data/B001_E.sac', 'data/B001_N.sac', 'data/B001_Z.sac')
        self.average_data = run_model(self.data, self.param, plot=False, verbose=False, write=False)

    def test_num_windows(self):
        self.assertEqual(self.average_data.num_windows, 100)

    def test_frequencies(self):
        self.setB = (3.3431085043988276, 0.8993157380254146, -3.3431085043988276, 16.53958944281525)
        frequencies = self.average_data.get_frequency(20)
        for truth, calculated in zip(self.setB, frequencies):
            self.assertAlmostEqual(truth, calculated)

    def test_plot(self):
        from hvarma import plot_hvratio
        self.setCa = [-20.,         -19.96089932, -19.92179863, -19.88269795, -19.84359726,
                      -19.80449658, -19.76539589, -19.72629521, -19.68719453, -19.64809384]
        self.setCb = [0.81721816, 0.8062195,  0.78768426, 0.77677296, 0.77333868, 0.76196565,
                      0.75378387, 0.72431371, 0.70199212, 0.67628692]

        fig = plot_hvratio(self.average_data, self.param, write_png=False)
        plt_data = fig.axes[0].get_lines()[0].get_data()
        assert_array_almost_equal(plt_data[0][:10], self.setCa)
        assert_array_almost_equal(plt_data[1][:10], self.setCb)

    def test_hvarma_call(self):
        from hvarma import HVarma as HVarma_test
        from hvarma.processing import HVarma
        self.assertIs(HVarma_test, HVarma)


if __name__ == '__main__':
    unittest.main(verbosity=2)