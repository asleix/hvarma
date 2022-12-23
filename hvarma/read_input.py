from dataclasses import dataclass, fields, field
from collections.abc import Sequence
import numpy as np
import os


@dataclass
class ArmaParam:
    """ Parameters for computation """
    model_order: int = 30
    maxtau:      int = 128
    mu:          float = 0.5
    nu:          float = 0.5
    nfir:        int = 40
    neg_freq:    float = -20
    pos_freq:    float = 20
    freq_points: int = 1024
    window_size: int = 512
    overlap:     int = 256
    max_windows: int = 1000
    freq_conf:   float = 20
    plot_conf:   float = 50
    output_dir:  str = '.'

    def __post_init__(self):
        """ Match the input values type during initialization. """
        for arg in fields(self):
            value = getattr(self, arg.name)
            if not isinstance(value, arg.type):
                setattr(self, arg.name, arg.type(value))

        self.assert_parameters()

    def assert_parameters(self):
        """ Perform some assertions to ensure parameters are consistent. """
        assert 2*self.maxtau <= self.window_size
        assert self.model_order <= self.maxtau
        assert self.overlap <= self.window_size
        assert self.neg_freq < self.pos_freq
        assert self.nfir <= self.maxtau

    @classmethod
    def from_dict(cls, args_dict):
        assert type(args_dict) is dict, \
            f"args_dict is not a dict but {type(args_dict)}"

        defaults = ArmaParam.get_fields_list()
        for arg in args_dict:
            if arg not in defaults:
                raise AttributeError(f'Parameter {arg} is unknown.')

        return cls(**args_dict)

    @classmethod
    def from_file(cls, args_file):
        assert os.path.exists(args_file), "Arguments file does not exist"
        default_values = ArmaParam.read_params(args_file)
        return cls.from_dict(default_values)

    @classmethod
    def get_fields_list(cls):
        return [arg.name for arg in fields(cls)]

    @staticmethod
    def read_params(filename):
        """ Read arguments from file. """
        values = {}
        with open(filename, 'r') as file:
            for line in file:
                try:
                    parameter, value = line.rstrip().split('=')
                except ValueError:
                    raise SyntaxError(f'Wrong syntax in {filename}. '
                                      f'Use: parameter=value')
                values[parameter] = value

        return values

    def get_dict(self):
        from dataclasses import asdict
        return asdict(self)

    def __copy__(self):
        return ArmaParam(self.get_dict())

    def copy(self):
        from copy import copy
        return copy(self)

    def update(self, args):
        assert type(args) is dict, "Usage: arg={'param':value}"
        dict_values = self.get_dict()

        for key, value in args.items():
            dict_values[key] = value

        return ArmaParam.from_dict(dict_values)


@dataclass
class Data:
    dataZ: Sequence
    dataN: Sequence
    dataE: Sequence
    sampling_rate: float
    station: str
    copy_data: bool = field(default=True, repr=False)

    def __post_init__(self):
        """ Match the input values type during initialization. """
        self.dataZ = np.array(self.dataZ, dtype=np.float64, copy=self.copy_data)
        self.dataN = np.array(self.dataN, dtype=np.float64, copy=self.copy_data)
        self.dataE = np.array(self.dataE, dtype=np.float64, copy=self.copy_data)
        self.sampling_rate = float(self.sampling_rate)
        self.station = str(self.station)

        assert len(self.dataE) == len(self.dataN) == len(self.dataZ), \
            "Data do not have the same size in Z, N or E directions"

        self.size = len(self.dataZ)

    @classmethod
    def from_sac(cls, Z_fname, N_fname, E_fname):
        """ Read data from SAC """
        from obspy import read

        def get_data(filename):
            st = read(filename, debug_headers=True)
            data = st[0].data
            return np.array(data, dtype=np.float64), st[0].stats

        dataE, headerE = get_data(E_fname)
        dataN, headerN = get_data(N_fname)
        dataZ, headerZ = get_data(Z_fname)

        attributes = ['sampling_rate', 'station']
        for key in attributes:
            if headerE[key] != headerN[key] or headerN[key] != headerZ[key]:
                raise ValueError('Fields have different header value at: ' + key)

        return cls(dataZ, dataN, dataE, headerE['sampling_rate'], headerE['station'])

    def make_window(self, start, size, copy=False):
        """ Create a slice of data. """
        assert start + size < self.size, "Window exceeds data size"
        dataZ = self.dataZ[start:start + size]
        dataN = self.dataN[start:start + size]
        dataE = self.dataE[start:start + size]
        return Data(dataZ, dataN, dataE, self.sampling_rate, self.station, copy_data=copy)

    def copy(self):
        # To be consistent with an old bug, this function returns itself!
        return self  # Data(self.dataZ, self.dataN, self.dataE, self.sampling_rate, self.station)
