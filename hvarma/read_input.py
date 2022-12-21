import numpy as np
from obspy import read
import os

DEFAULT_PARAMS = {
    'model_order':       int,
    'maxtau':            int,
    'mu':                float,
    'nu':                float,
    'nfir':              int,
    'neg_freq':          float,
    'pos_freq':          float,
    'freq_points':       int,
    'window_size':       int,
    'overlap':           int,
    'max_windows':       int,
    'freq_conf':         float,
    'plot_conf':         float,
    'output_dir':        str,
}


class ArmaParam:
    """ Parameters for computation """
    model_order: int
    maxtau:      int
    mu:          float
    nu:          float
    nfir:        int
    neg_freq:    float
    pos_freq:    float
    freq_points: int
    window_size: int
    overlap:     int
    max_windows: int
    freq_conf:   float
    plot_conf:   float
    output_dir:  str

    def __init__(self, arg='default'):
        assert isinstance(arg, str) or isinstance(arg, dict), "Invalid argument"
        filename = arg if isinstance(arg, str) else 'default'
        default_filename = os.path.dirname(os.path.abspath(__file__)) + \
                           os.path.sep + 'default_args.txt'
        default_values = self.read_params(default_filename)

        if filename == 'default':
            filename = default_filename
        self.filename = filename

        # Overwrite default values if provided in argument
        values = self.read_params(filename)
        default_values.update(values)
        if isinstance(arg, dict):
            default_values.update(arg)
        values = default_values

        # Set class attributes from dictionary
        for param, ptype in DEFAULT_PARAMS.items():
            values[param] = ptype(values[param])
            setattr(self, param, values[param])

        self.dict_values = values
        self.assert_parameters()

    @staticmethod
    def get_param_list():
        return list(DEFAULT_PARAMS.keys())

    def read_params(self, filename):
        """ Read arguments from file. """
        values = {}
        parameter_list = ArmaParam.get_param_list()
        with open(filename, 'r') as f:
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

        return values

    def assert_parameters(self):
        """ Perform some assertions to ensure parameters are consistent. """
        assert 2*self.maxtau <= self.window_size
        assert self.model_order <= self.maxtau
        assert self.overlap <= self.window_size
        assert self.neg_freq < self.pos_freq
        assert self.nfir <= self.maxtau

    def get_params(self):
        return self.dict_values.copy()

    def __copy__(self):
        return ArmaParam(self.get_params())

    def copy(self):
        from copy import copy
        return copy(self)

    def update(self, arg):
        assert type(arg) is dict, "Usage: arg={'param':value}"
        dict_values = self.get_params()
        assert all([key in dict_values for key in arg]), "Wrong parameter"

        for key, value in arg.items():
            dict_values[key] = value

        return ArmaParam(dict_values)


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


