import sys
import os
import time
from collections import OrderedDict
import numpy as np
from hvarma.read_input import Window
from hvarma.processing import HVarma, AverageData
from hvarma.write_output import progress_bar, write_results, plot_hvratio, plot_order_search


def get_data_windows(data, size, overlap):
    """ Generator of data slices from data of a given size,
        overlapping one another """
    if data.size < size:
        raise Exception('Window exceeds available data')

    for start in range(0, data.size-size+1, size-overlap):
        yield Window(data, start, size)


def run_model(data, param, plot=False, verbose=True, write=False):
    out = sys.stdout if verbose else open(os.devnull, "w")

    # Start windowing
    beg = time.time()
    processed_windows = []

    progress = progress_bar(data.size, param.window_size, param.overlap, param.max_windows, file=out)
    for idx, data_window in enumerate(get_data_windows(data, param.window_size, param.overlap)):
        next(progress)
        model = HVarma(data_window, param)
        model.solve_arma()
        model.get_coherence()

        processed_windows.append(model)
        if idx + 1 == param.max_windows:  # Limited windows version
            break

    print('Elapsed:', round((time.time() - beg) / 60, 1), 'min', file=out)

    print('Retrieving spectra...', file=out)
    results = AverageData(processed_windows, param)

    if write:
        print('Constructing output...', file=out)
        write_results(data, param, results)

    if plot:
        plot_hvratio(param, results, write=True)

    pos_freq, pos_err, neg_freq, neg_err = results.get_frequency(param.freq_conf)
    print('Estimated positive resonance frequency: {:.6f} Hz, error: {:.6f} Hz'.format(pos_freq, pos_err), file=out)
    print('Estimated negative resonance frequency: {:.6f} Hz, error: {:.6f} Hz'.format(neg_freq, neg_err), file=out)

    if not verbose:
        out.close()

    return results


def get_difference(cur, prev):
    """ Subtract current and previous frequencies (positive, negative) """
    freqs_cur = cur.get_frequency(20)
    freqs_prev = prev.get_frequency(20)

    pos_diff = freqs_cur[0]-freqs_prev[0]
    neg_diff = freqs_cur[2]-freqs_prev[2]

    return pos_diff, neg_diff


def convergence_condition(pos_diff, neg_diff, tol):
    """ Condition to determine if a given model order has converged within tolerance """
    return abs(pos_diff)+abs(neg_diff) < 2*tol


def get_results_for_order(data, param, tested_orders, order):
    """ Run model for given order and order-3
        if not already computed in tested_orders."""
    param_cur = param.update({'model_order': order})
    param_prev = param.update({'model_order': order-3})

    if order not in tested_orders:
        tested_orders[order] = run_model(data, param_cur, plot=False, verbose=False, write=False)

    if order-3 not in tested_orders:
        tested_orders[order-3] = run_model(data, param_prev, plot=False, verbose=False, write=False)

    return tested_orders[order], tested_orders[order-3]


def is_converged(data, param, tested_orders, order, tol=0.1):
    """ Check if a model order is sufficient to model given data (convergence criterion) """
    results_cur, results_prev = get_results_for_order(data, param, tested_orders, order)
    pos_diff, neg_diff = get_difference(results_cur, results_prev)
    converged = convergence_condition(pos_diff, neg_diff, tol=tol)
    return converged


def binary_search(data, param, tested_orders, low_p, high_p, tol=0.1, verbose=False):
    """ Find smallest converged order in range low_p, high_p """
    out = sys.stdout if verbose else open(os.devnull, "w")
    print('Refining order within found bounds:', end='', file=out)
    while low_p < high_p:

        mid_p = (high_p + low_p) // 2
        print(f' {mid_p}', end='', file=out)
        sys.stdout.flush()
        if is_converged(data, param, tested_orders, mid_p, tol):
            high_p = mid_p
        else:
            low_p = mid_p+1

    print(file=out)
    if not verbose:
        out.close()
    assert low_p == high_p, f'{low_p} != {high_p}'
    return low_p


def find_optimal_order_fast(data, param, tol=0.05, start_order=4, output_dir='.',
                            plot=False, verbose=False, write=False):
    """
    Use fast algorithm to find a small converged hvarma order for given data.
    """
    out = sys.stdout if verbose else open(os.devnull, "w")
    assert start_order >= 4

    beg = time.time()
    order = start_order
    tested_orders = OrderedDict()

    print('Finding order upper bound. Tested orders:', end='', file=out)
    sys.stdout.flush()

    while not is_converged(data, param, tested_orders, order, tol=tol):
        print(f' {order}', end='', file=out)
        sys.stdout.flush()
        if order == param.maxtau:
            break
        order = int(order * 2)
        if order > param.maxtau:
            order = param.maxtau

    print(file=out)
    # Now bisection search to refine order
    final_order = binary_search(data, param, tested_orders, int(order / 2), order,
                                tol=tol, verbose=verbose)
    if plot:
        plot_order_search(tested_orders, final_order, data.station, tol=tol, output_dir=output_dir)

    print('Elapsed:', round((time.time() - beg) / 60, 1), 'min', file=out)

    result = is_converged(data, param, tested_orders, final_order, tol=tol)
    if result:
        print('Final order', final_order, file=out)

    if not verbose:
        out.close()
    if not result:
        raise RuntimeError('Could not find convergence with the given parameters')

    return final_order


def find_optimal_order(data, param, tol=0.05, start_order=4, output_dir='.',
                       plot=False, verbose=False, write=False, method='fast'):
    """
    Find a small hvarma order that suffices to describe data.
    The returned order satisfies a convergence criterion.
    """
    if method == 'fast':
        return find_optimal_order_fast(data, param, tol=tol, start_order=start_order,
                                       output_dir=output_dir,
                                       plot=plot, verbose=verbose, write=write)

    assert 0, f"Method {method} not available"
