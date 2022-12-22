import argparse
from hvarma import Data, ArmaParam, run_model, write_results, plot_hvratio


def select_parameters_from_args(args):
    """ Filter out ArmaParam attributes from commandline arguments """
    args_dict = {}
    for arg in vars(args):
        if arg in ArmaParam.get_param_list():
            if getattr(args, arg) is not None:
                args_dict[arg] = getattr(args, arg)
    return args_dict


def write_frequencies_in_file(param, results):
    """ Write the frequency results in a file """
    from hvarma.write_output import generate_filename
    filename = generate_filename(param.output_dir, results.station,
                                 param.model_order, results.num_windows)
    pos_freq, pos_err, neg_freq, neg_err = results.get_frequency(param.freq_conf)
    with open(filename + '.res', 'w') as fres:
        fres.write('POSITIVE {:.6f} error {:.6f} Hz '
                   'NEGATIVE {:.6f} error {:.6f} Hz\n'
                   .format(pos_freq, pos_err, neg_freq, neg_err))


def main(args):
    data = Data(Z_fname=args.Z_fname,
                N_fname=args.N_fname,
                E_fname=args.E_fname)

    param = ArmaParam(args.args_file)
    args_dict = select_parameters_from_args(args)
    param = param.update(args_dict)

    if not args.silent:
        print('Data read correctly')

    results = run_model(data, param, verbose=not args.silent)
    write_results(data, param, results)
    plot_hvratio(param, results, format='png')
    write_frequencies_in_file(param, results)


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
    parser.add_argument('--args_file', type=str, help="File with all the default arguments.", default='default')
    parser.add_argument('--freq_points', type=int, help="Number of frequency points to "
                                                        "calculate between neg_freq and pos_freq")
    parser.add_argument('--freq_conf', type=float, help='Frequency confidence interval.')
    parser.add_argument('--silent', help="No output to stdout", action='store_false', default=False)

    main(parser.parse_args())
