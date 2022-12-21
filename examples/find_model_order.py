import argparse
from hvarma import Data, ArmaParam, find_optimal_order


def select_parameters_from_args(args):
    args_dict = {}
    for arg in vars(args):
        if arg in ArmaParam.get_param_list():
            if getattr(args, arg) is not None:
                args_dict[arg] = getattr(args, arg)
    return args_dict


def main(args):

    data = Data(Z_fname=args.Z_fname,
                N_fname=args.E_fname,
                E_fname=args.N_fname)
    param = ArmaParam(args.args).update({
        'freq_points': 2000,
        'neg_freq':  -10,
        'pos_freq':  10,
        'max_windows': 250,
        'window_size': 512,
    })
    args_dict = select_parameters_from_args(args)
    param = param.update(args_dict)
    print(param.get_params())
    print('Data read correctly')
    find_optimal_order(data, param, 0.05, start_order=4, output_dir='.', verbose=True, plot=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('Z_fname', type=str, help="SAC data in direction Z")
    parser.add_argument('N_fname', type=str, help="SAC data in direction N")
    parser.add_argument('E_fname', type=str, help="SAC data in direction E")
    parser.add_argument('--model_order', type=int, help="HVARMA model order")
    parser.add_argument('--output_dir', type=str, help="Directory in which to store output")
    parser.add_argument('--max_windows', type=int, help="Maximum number of windows to explore within data.")
    parser.add_argument('--window_size', type=int, help="Size of each individual window.")
    parser.add_argument('--overlap', type=int, help="Overlap between individual windows.")
    parser.add_argument('--args', type=str, help="File with all the default arguments.", default='default')
    parser.add_argument('--freq_points', type=int, help="Number of frequency points to "
                                                        "calculate between neg_freq and pos_freq")
    parser.add_argument('--silent', help="No output to stdout", action='store_false', default=False)
    main(parser.parse_args())
