import unittest
import io
import contextlib
import os, sys
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))


from include.read_input import ArmaParam, Data, Window
from include.write_output import progress_bar, write_data, pretty_show


class DataTest(unittest.TestCase):

    def setUp(self):
        self.Z_fn = 'data/B001_Z.sac'
        self.N_fn = 'data/B001_N.sac'
        self.E_fn = 'data/B001_E.sac'

    def test_input1(self):
        data = Data(Z_fname=self.Z_fn, N_fname=self.N_fn, E_fname=self.E_fn)

        self.assertEqual(data.size, 1886784)
        self.assertAlmostEqual(data.sampling_rate, 100.0, places=8)
        self.assertEqual(str(data.starttime), '2019-03-06T22:23:33.103400Z')
        self.assertEqual(str(data.endtime), '2019-03-07T03:38:00.933400Z')

    def test_wrong(self):
        with self.assertRaises(Exception):
            data = Data(Z_fname=self.Z_fn)

        with self.assertRaises(Exception):
            data = Data(E_fname=self.E_fn)

        with self.assertRaises(Exception):
            data = Data(F_fname=self.F_fn)

    def test_window(self):
        data = Data(Z_fname=self.Z_fn, N_fname=self.N_fn, E_fname=self.E_fn)
        win = Window(data, 0, 100)
        from numpy.testing import assert_array_almost_equal
        assert_array_almost_equal(win.dataE, data.dataE[0:100])

        with self.assertRaises(AttributeError):
            Window('wrong_arg', 0, 10)


class ArmaParamTest(unittest.TestCase):

    def setUp(self):
        self.filename = 'test/resources/args1.txt'

    def test_input1(self):
        param = ArmaParam(filename=self.filename)

        self.assertEqual(param.p, 74)
        self.assertEqual(param.maxtau, 128)
        self.assertEqual(param.mu, 0.5)
        self.assertEqual(param.nu, 0.5)
        self.assertEqual(param.nfir, 40)
        self.assertEqual(param.f0, -20)
        self.assertEqual(param.f1, 20)
        self.assertEqual(param.npun, 1024)
        self.assertEqual(param.wsize, 512)
        self.assertEqual(param.overlap, 256)
        self.assertEqual(param.maxwin, 1000)
        self.assertEqual(param.freq_conf, 20)
        self.assertEqual(param.plot_conf, 50)
        self.assertEqual(param.oname, 'default')


class ProgressBarTest(unittest.TestCase):
    def test_bar_multi(self):
        f = io.StringIO()
        with contextlib.redirect_stdout(f): # catch stdout
            # progress_bar(size, window_size, overlap, max_win)

            progress = progress_bar(100, 10, 1, 20)
            count = sum([1 for idx in progress])
            self.assertEqual(count, 11) # Simple Overlap

            progress = progress_bar(100, 10, 0, 20)
            count = sum([1 for idx in progress])
            self.assertEqual(count, 10) # No overlap

            progress = progress_bar(90, 10, 5, 30)
            count = sum([1 for idx in progress])
            self.assertEqual(count, 17) # More complex overlap

            with self.assertRaises(AssertionError): # Wrong parameters
                progress = progress_bar(200, 10, 50, 110)
                count = sum([1 for idx in progress])

            progress = progress_bar(200, 10, 1, 4)
            count = sum([1 for idx in progress])
            self.assertEqual(count, 4) # max windows parameter


if __name__ == '__main__':
    unittest.main(verbosity=2)

