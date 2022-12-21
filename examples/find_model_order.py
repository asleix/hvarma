import argparse
from hvarma import Data, ArmaParam, find_optimal_order


def main(args):

    data = Data(Z_fname=args.Z_fname,
                N_fname=args.E_fname,
                E_fname=args.N_fname)
    params = ArmaParam({
            'arma_order': 10,
            'ini_freq': -10,
            'fin_freq': 10,
            'num_points': 2000,
            'window_size': 1024,
            'max_windows': 50
    })

    print('Data read correctly')
    find_optimal_order(data, params, 0.1, start_order=4, output_dir='.', verbose=True, plot=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('Z_fname', type=str)
    parser.add_argument('N_fname', type=str)
    parser.add_argument('E_fname', type=str)

    main(parser.parse_args())
