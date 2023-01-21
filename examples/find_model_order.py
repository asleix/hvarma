"""
Copyright (c) 2022, Spanish National Research Council (CSIC)

Script to run HVarma of a specific order on input data.

HVarma estimates the transfer function in a surface layer for
three-dimensional micro-tremor seismogram data.

Usage example:
    python run.py Z_data.sac N_data.sac E_data.sac --start_order=10
"""

import argparse
from hvarma import Data, ArmaParam, find_optimal_order, plot_order_search


def select_parameters_from_args(args):
    args_dict = {}
    for arg in vars(args):
        if arg in ArmaParam.get_fields_list():
            if getattr(args, arg) is not None:
                args_dict[arg] = getattr(args, arg)
    return args_dict


def main(args):
    data = Data.from_sac(Z_fname=args.Z_fname,
                         N_fname=args.E_fname,
                         E_fname=args.N_fname)
    if args.args is not None:
        param = ArmaParam.from_file(args.args)
    else:
        param = ArmaParam.from_dict({
            'freq_points': 2000,
            'maxtau':  64,
            'neg_freq':    -10,
            'pos_freq':    10,
            'max_windows': 10,
            'window_size': 512,
        })
    args_dict = select_parameters_from_args(args)
    param = param.update(args_dict)
    if not args.silent:
        print('Data read correctly')
    results = find_optimal_order(data, param, 0.05, start_order=args.start_order, verbose=not args.silent)
    plot_order_search(results, param.output_dir)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('Z_fname', type=str, help="SAC data in direction Z")
    parser.add_argument('N_fname', type=str, help="SAC data in direction N")
    parser.add_argument('E_fname', type=str, help="SAC data in direction E")
    parser.add_argument('--start_order', type=int, help="HVARMA model order", default=4)
    parser.add_argument('--output_dir', type=str, help="Directory in which to store output")
    parser.add_argument('--max_windows', type=int, help="Maximum number of windows to explore within data.")
    parser.add_argument('--window_size', type=int, help="Size of each individual window.")
    parser.add_argument('--overlap', type=int, help="Overlap between individual windows.")
    parser.add_argument('--args', type=str, help="File with all the default arguments.", default=None)
    parser.add_argument('--freq_points', type=int, help="Number of frequency points to "
                                                        "calculate between neg_freq and pos_freq")
    parser.add_argument('--silent', help="No output to stdout", action='store_false', default=False)
    main(parser.parse_args())
