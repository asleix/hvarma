import time
import numpy as np
from hvarma.read_input import ArmaParam, Data, Window
from hvarma.processing import HVarma, AverageData
from hvarma.write_output import progress_bar, write_data, pretty_show


def get_data_windows(data, size, overlap):
    """ Generator of data slices from data of a given size,
        overlapping one another """
    if data.size < size:
        raise Exception('Window exceeds available data')

    for start in range(0, data.size-size+1, size-overlap):
        yield Window(data, start, size)


def main(args):
    # Read data and parameters
    param = ArmaParam()
    data = Data(Z_fname=args.Z_fname, N_fname=args.E_fname, E_fname=args.N_fname)
    print('Data read correctly')

    # Start windowing
    beg = time.time()
    processed_windows = []
    nwin = param.maxwin
    progress = progress_bar(data.size, param.wsize, param.overlap, param.maxwin)
    for idx, data_window in enumerate(get_data_windows(data, param.wsize, param.overlap)):
        nwin = next(progress)
        model = HVarma(data_window, param)
        model.solve_arma()
        model.get_coherence()

        processed_windows.append(model)
        if idx+1 == param.maxwin:  # Limited windows version
            break

    print('Elapsed:', round((time.time()-beg)/60, 1), 'min')

    print('Retrieving spectra...')
    compute_averages = AverageData(processed_windows, param)

    print('Constructing output...')
    if param.oname == 'default':
        outfile = 'output/{}_p{}_win{}'.format(data.station, param.p, nwin)
    else:
        if not param.oname.endswith('/'):
            param.oname += '/'
        outfile = param.oname + '{}_p{}_win{}'.format(data.station, param.p, nwin)
    err = param.plot_conf/2
    write_data(np.linspace(param.f0, param.f1, param.npun), compute_averages.get_spectrum_percentile(50),
               compute_averages.get_spectrum_percentile(50-err), compute_averages.get_spectrum_percentile(50+err),
               compute_averages.get_coherence_percentile(50), outfile+'.txt')

    pos_freq, pos_err, neg_freq, neg_err = compute_averages.get_frequency(param.freq_conf)
    pretty_show(param.f0, param.f1, param.npun, compute_averages.get_spectrum_percentile(50),
                compute_averages.get_spectrum_percentile(50-err), compute_averages.get_spectrum_percentile(50+err),
                compute_averages.get_coherence_percentile(50), data.station,
                pos_freq, pos_err, neg_freq, neg_err, filename=outfile+'.png')

    print('Estimated positive resonance frequency: {:.6f} Hz, error: {:.6f} Hz'.format(pos_freq, pos_err))
    print('Estimated negative resonance frequency: {:.6f} Hz, error: {:.6f} Hz'.format(neg_freq, neg_err))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('Z_fname', type=str)
    parser.add_argument('N_fname', type=str)
    parser.add_argument('E_fname', type=str)

    main(parser.parse_args())
