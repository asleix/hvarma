import unittest
import io
import contextlib


class DataTest(unittest.TestCase):

    def setUp(self):
        self.Z_fn = 'test/resources/B001_Z.sac'
        self.N_fn = 'test/resources/B001_N.sac'
        self.E_fn = 'test/resources/B001_E.sac'

    def test_input1(self):
        from hvarma.read_input import Data
        data = Data.from_sac(Z_fname=self.Z_fn, N_fname=self.N_fn, E_fname=self.E_fn)

        self.assertEqual(data.size, 1886784)
        self.assertAlmostEqual(data.sampling_rate, 100.0, places=8)

    def test_wrong(self):
        from hvarma.read_input import Data
        with self.assertRaises(Exception) as msg:
            data = Data.from_sac(Z_fname=self.Z_fn)
        self.assertIn('positional arguments', str(msg.exception))

        with self.assertRaises(Exception) as msg:
            data = Data.from_sac(E_fname=self.E_fn)
        self.assertIn('positional arguments', str(msg.exception))

        with self.assertRaises(Exception) as msg:
            data = Data.from_sac(N_fname=self.N_fn)
        self.assertIn('positional arguments', str(msg.exception))

    def test_window(self):
        from hvarma.read_input import Data
        data = Data.from_sac(Z_fname=self.Z_fn, N_fname=self.N_fn, E_fname=self.E_fn)
        win = data.make_window(0, 100)
        from numpy.testing import assert_array_almost_equal
        assert_array_almost_equal(win.dataE, data.dataE[0:100])


class ArmaParamTest(unittest.TestCase):

    def setUp(self):
        self.filename = 'test/resources/args1.txt'
        self.filename2 = 'test/resources/args2.txt'

    def test_input1(self):
        from hvarma.read_input import ArmaParam
        param = ArmaParam.from_file(self.filename)

        self.assertEqual(param.model_order, 74)
        self.assertEqual(param.maxtau, 128)
        self.assertEqual(param.mu, 0.5)
        self.assertEqual(param.nu, 0.5)
        self.assertEqual(param.nfir, 40)
        self.assertEqual(param.neg_freq, -20)
        self.assertEqual(param.pos_freq, 20)
        self.assertEqual(param.freq_points, 1024)
        self.assertEqual(param.window_size, 512)
        self.assertEqual(param.overlap, 256)
        self.assertEqual(param.max_windows, 100)
        self.assertEqual(param.freq_conf, 20)
        self.assertEqual(param.plot_conf, 50)
        self.assertEqual(param.output_dir, '.')

    def test_update(self):
        from hvarma.read_input import ArmaParam

        self.setA = ArmaParam.from_file(self.filename2)
        param = ArmaParam.from_file(self.filename).update({
            'model_order': 10,
            'neg_freq': -10,
            'pos_freq': 10,
            'freq_points': 2000,
            'window_size': 512,
            'max_windows': 250
        })
        self.assertDictEqual(self.setA.get_dict(), param.get_dict())


class ProgressBarTest(unittest.TestCase):
    def test_bar_multi(self):
        from hvarma.write_output import progress_bar
        f = io.StringIO()
        with contextlib.redirect_stdout(f):  # catch stdout
            # progress_bar(size, window_size, overlap, max_win)

            progress = progress_bar(100, 10, 1, 20)
            count = sum([1 for _ in progress])
            self.assertEqual(count, 11)  # Simple Overlap

            progress = progress_bar(100, 10, 0, 20)
            count = sum([1 for _ in progress])
            self.assertEqual(count, 10)  # No overlap

            progress = progress_bar(90, 10, 5, 30)
            count = sum([1 for _ in progress])
            self.assertEqual(count, 17)  # More complex overlap

            with self.assertRaises(AssertionError):  # Wrong parameters
                progress = progress_bar(200, 10, 50, 110)
                sum([1 for _ in progress])

            progress = progress_bar(200, 10, 1, 4)
            count = sum([1 for _ in progress])
            self.assertEqual(count, 4)  # max windows parameter


if __name__ == '__main__':
    unittest.main(verbosity=2)
