import numpy as np
from obspy import read
import os


class ArmaParam:
    """ Parameters for computation """
    def __init__(self, filename='default'):
        if filename == 'default':
            filename = os.path.dirname(os.path.abspath(__file__)) + os.path.sep + '../args.txt'
        self.filename = filename
        values = self.read_params()

        self.p = int(values['arma_order'])
        self.maxtau = int(values['maxtau'])
        self.mu = float(values['forward_weight'])
        self.nu = float(values['backward_weight'])
        self.nfir = int(values['nfir'])
        self.f0 = float(values['ini_freq'])
        self.f1 = float(values['fin_freq'])
        self.npun = int(values['num_points'])
        self.wsize = int(values['window_size'])
        self.overlap = int(values['window_shift'])
        self.maxwin = int(values['max_windows'])
        self.freq_conf = float(values['freq_conf_interval'])
        self.plot_conf = float(values['spectral_conf_interval'])
        self.oname = str(values['output_path'])

        self.assert_parameters()

    def read_params(self):
        """ Read arguments from file. """
        values = {}
        parameter_list = ['arma_order', 'maxtau', 'nfir', 'ini_freq', 'fin_freq',
                          'num_points', 'window_size', 'window_shift', 'max_windows',
                          'freq_conf_interval', 'spectral_conf_interval', 'output_path',
                          'forward_weight', 'backward_weight']
        with open(self.filename, 'r') as f:
            for line in f:
                line = line.replace(' ', '')
                line = line.replace('\n', '')
                try:
                    parameter, value = line.split('=')
                except ValueError:
                    raise SyntaxError('Wrong syntax. Use: parameter=value')
                if parameter not in parameter_list:
                    raise AttributeError('Parameter ' + parameter + ' unknown.')
                values[parameter] = value

        if len(values.keys()) != len(parameter_list):
            raise AttributeError('Wrong number of input parameters. Check arguments file.')

        return values

    def assert_parameters(self):
        """ Perform some assertions to ensure parameters are consistent. """
        assert 2*self.maxtau <= self.wsize
        assert self.p <= self.maxtau
        assert self.overlap <= self.wsize
        assert self.f0 < self.f1
        assert self.nfir <= self.maxtau


class Data:
    """ Main class for data storage """
    def __init__(self, E_fname=None, N_fname=None, Z_fname=None):
        """ Initialize class variables """
        self.dataE = None
        self.dataN = None
        self.dataZ = None
        self.size = None
        self.sampling_rate = None
        self.starttime = None
        self.endtime = None
        self.station = None

        # SAC initialization
        assertion = [E_fname is not None, N_fname is not None, Z_fname is not None]
        if all(assertion):
            self.read_sac(E_fname, N_fname, Z_fname)
        elif any(assertion):
            raise Exception('Bad data initialization, missing arguments')

    def read_sac(self, E_fname, N_fname, Z_fname):
        """ Read data from SAC """
        def get_data(filename):
            st = read(filename, debug_headers=True)
            data = st[0].data
            return np.array(data, dtype=np.float64), st[0].stats

        self.dataE, headerE = get_data(E_fname)
        self.dataN, headerN = get_data(N_fname)
        self.dataZ, headerZ = get_data(Z_fname)

        def is_data_sync(headerE, headerN, headerZ):
            fields = ['sampling_rate', 'starttime', 'endtime']
            for field in fields:
                if headerE[field] != headerN[field] or headerN[field] != headerZ[field]:
                    raise Exception('Data not sync at ' + field)

        is_data_sync(headerE, headerN, headerZ)

        self.sampling_rate = headerE['sampling_rate']
        self.starttime = headerE['starttime']
        self.endtime = headerE['endtime']
        self.station = headerE['station']

        if len(self.dataE) != len(self.dataN) or len(self.dataN) != len(self.dataZ):
            raise Exception('Input files have different sizes')

        self.size = len(self.dataE)


class Window(Data):
    """ Subclass that slices a Window from Data """
    def __init__(self, data, start, size):
        """ data is Data object, size integer, size integer """
        # split data window
        super().__init__()
        if not isinstance(data, Data):
            raise AttributeError('Wrong "data" argument.')

        self.dataE = data.dataE[start:start+size]
        self.dataN = data.dataN[start:start+size]
        self.dataZ = data.dataZ[start:start+size]
        self.size = size
        self.sampling_rate = data.sampling_rate
        self.starttime = data.starttime
        self.endtime = data.endtime
        self.station = data.station


