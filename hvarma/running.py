import sys
import os
import time
import numpy as np
from hvarma.read_input import Window
from hvarma.processing import HVarma, AverageData
from hvarma.write_output import progress_bar, write_data, pretty_show


def get_data_windows(data, size, overlap):
    """ Generator of data slices from data of a given size,
        overlapping one another """
    if data.size < size:
        raise Exception('Window exceeds available data')

    for start in range(0, data.size-size+1, size-overlap):
        yield Window(data, start, size)


def run_model(data, param, plot=False, verbose=False, write=False):
    out = sys.stdout if verbose else open(os.devnull, "w")

    # Start windowing
    beg = time.time()
    processed_windows = []

    progress = progress_bar(data.size, param.wsize, param.overlap, param.maxwin, file=out)
    for idx, data_window in enumerate(get_data_windows(data, param.wsize, param.overlap)):
        next(progress)
        model = HVarma(data_window, param)
        model.solve_arma()
        model.get_coherence()

        processed_windows.append(model)
        if idx + 1 == param.maxwin:  # Limited windows version
            break

    print('Elapsed:', round((time.time() - beg) / 60, 1), 'min', file=out)

    print('Retrieving spectra...', file=out)
    compute_averages = AverageData(processed_windows, param)

    print('Constructing output...', file=out)
    num_windows = len(processed_windows)
    if param.oname == 'default':
        outfile = 'output/{}_p{}_win{}'.format(data.station, param.p, num_windows)
    else:
        if not param.oname.endswith('/'):
            param.oname += '/'
        outfile = param.oname + '{}_p{}_win{}'.format(data.station, param.p, num_windows)
    err = param.plot_conf / 2

    if write:
        write_data(np.linspace(param.f0, param.f1, param.npun), compute_averages.get_spectrum_percentile(50),
                   compute_averages.get_spectrum_percentile(50 - err),
                   compute_averages.get_spectrum_percentile(50 + err),
                   compute_averages.get_coherence_percentile(50), outfile + '.txt')

    pos_freq, pos_err, neg_freq, neg_err = compute_averages.get_frequency(param.freq_conf)

    if plot:
        pretty_show(param.f0, param.f1, param.npun, compute_averages.get_spectrum_percentile(50),
                    compute_averages.get_spectrum_percentile(50 - err),
                    compute_averages.get_spectrum_percentile(50 + err),
                    compute_averages.get_coherence_percentile(50), data.station,
                    pos_freq, pos_err, neg_freq, neg_err, filename=outfile + '.png')

    print('Estimated positive resonance frequency: {:.6f} Hz, error: {:.6f} Hz'.format(pos_freq, pos_err), file=out)
    print('Estimated negative resonance frequency: {:.6f} Hz, error: {:.6f} Hz'.format(neg_freq, neg_err), file=out)

    if not verbose:
        out.close()

    return compute_averages
