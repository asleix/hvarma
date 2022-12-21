import argparse
from hvarma import Data, ArmaParam, run_model


def main(args):
    # Use default parameters
    param = ArmaParam()
    data = Data(Z_fname=args.Z_fname,
                N_fname=args.E_fname,
                E_fname=args.N_fname)
    print('Data read correctly')
    run_model(data, param, plot=True, verbose=True, write=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('Z_fname', type=str)
    parser.add_argument('N_fname', type=str)
    parser.add_argument('E_fname', type=str)

    main(parser.parse_args())
