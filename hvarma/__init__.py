"""
Copyright (c) 2022, Spanish National Research Council (CSIC)

This source code is subject to the terms of the
GNU Lesser General Public License.

HVarma estimates the transfer function in a surface layer for
three-dimensional micro-tremor seismogram data.

Usage example:
    # assume data and param are Data and Param instances
    results = run_model(data, param, verbose=False)
    pos_freq, pos_err, neg_freq, neg_err = results.get_frequency(param.freq_conf)

    results = find_optimal_order(data, param, 0.05, start_order=4)
    print('Found order:', results.final_order)
"""

from .running import run_model, find_optimal_order
from .processing import HVarma, AverageData
from .read_input import Data, ArmaParam
from .write_output import plot_hvratio, write_results, plot_order_search
