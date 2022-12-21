import os, sys
import time
import numpy as np


def progress_bar(size, wsize, overlap, maxwin, file=sys.stdout):
    """ Return an instance of a generator to display the progress of window processing.
         Usage call: next(progress_bar_object).     Yields the current window number.   """
    assert wsize > overlap, "Window size must be larger than overlap"

    start = time.time()

    def show(j):  # Print or update progress bar on screen
        x = int(30 * j / nwin)
        now = time.time()-start
        time_left = np.ceil(now * (nwin-j)/j/60)
        file.write("%s[%s%s] %i/%i, time left: %i min\r" %
                   ('Progress: ', "#" * x, "." * (30 - x), j, nwin, time_left))
        file.flush()

    nwin = min((size - overlap) // (wsize - overlap), maxwin)

    for i in range(1, nwin):
        show(i)
        yield i

    show(nwin)  # Last iteration writes new line
    file.write("\n")
    file.flush()
    yield nwin


def pretty_show(f0, f1, npun, spectrum, low_err, upp_err, coherence, stat_name,
                pos_freq, pos_err, neg_freq, neg_err, filename='out.png'):
    """ Draw H/V ratio spectrum. It's a scatter plot H/V amplitude vs frequency.
        Display the spectrum along with upper and lower error bounds.
        Also display the estimated frequency of the peaks, both positive and negative. """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    freqs = np.linspace(f0, f1, npun)

    # Set Times font
    plt.rc('font', family='serif', size=12)
    plt.rcParams['mathtext.fontset'] = 'dejavuserif'

    # Create figure and axis
    fig, ax1 = plt.subplots()
    fig.set_size_inches(8,5)
    ax2 = ax1.twinx()
    ax1.set_xlabel('Frequency (Hz)', fontsize=14)
    ax1.set_xlim(f0, f1)
    ax1.set_ylabel('H/V', fontsize=14)
    ax2.set_ylabel('Coherence', fontsize=14)
    ax2.set_ylim(0, 1)

    # Plot spectrum and error bounds
    ax1.fill_between(freqs, low_err, upp_err, alpha=0.3, color='k')
    ax2.plot(freqs, coherence, '--', color='k', linewidth=1)
    ax1.plot(freqs, spectrum, 'o', mec='k', mfc='none', ms=4)

    # Annotate the station name and the peak frequencies
    ax2.text(0.85*f1 + 0.15*f0, 0.88, stat_name, fontweight='bold', fontsize=14)
    if f1 > 0:
        pos_text = '$f_+ = {:.2f} \pm {:.2f}$ Hz'.format(pos_freq, pos_err)
        ax2.text(0.95*f0 + 0.05*f1, 0.92, pos_text, fontsize=12)
    if f0 < 0:
        neg_text = '$f_- = {:.2f} \pm {:.2f}$ Hz'.format(neg_freq, neg_err)
        ax2.text(0.95*f0 + 0.05*f1, 0.84, neg_text, fontsize=12)

    # Save figure
    if filename is not None:
        fig.savefig(filename, dpi=300, bbox_inches='tight')
    return fig


def plot_hvratio(average_data, param, write_png=True):
    """ Display the results of an hvarma model calculation.
        Draw H/V ratio spectrum. It's a scatter plot H/V amplitude vs frequency.
        Display the spectrum along with upper and lower error bounds.
        Also display the estimated frequency of the peaks, both positive and negative.
    """
    station_name = average_data.station
    num_windows = average_data.num_windows
    if write_png:
        if param.output_dir == 'default':
            outfile = './{}_p{}_win{}'.format(station_name, param.model_order, num_windows)
        else:
            if not param.output_dir.endswith('/'):
                param.output_dir += '/'
            outfile = param.output_dir + '{}_p{}_win{}'.format(station_name, param.model_order, num_windows)
        outfile += '.png'
    else:
        outfile = None

    err = param.plot_conf/2
    pos_freq, pos_err, neg_freq, neg_err = average_data.get_frequency(param.freq_conf)

    return pretty_show(param.neg_freq, param.pos_freq, param.freq_points, average_data.get_spectrum_percentile(50),
                       average_data.get_spectrum_percentile(50-err),
                       average_data.get_spectrum_percentile(50+err),
                       average_data.get_coherence_percentile(50), station_name,
                       pos_freq, pos_err, neg_freq, neg_err, filename=outfile)


def write_data(freqs, spectrum, low_err, up_err, coherence, filename='out.txt'):
    """ Write processed data in a text file. """
    with open(filename, 'w') as f:
        print('Frequency', 'H/V', 'Low_err', 'Upp_err', 'Coherence', file=f, sep=' ')
        for freq, spec, low, up, coh in zip(freqs, spectrum, low_err, up_err, coherence):
            print(freq, spec, low, up, coh, file=f, sep=' ')


def plot_order_search(orders, found_p, stat_name, tol=0.1, output_dir='.'):
    """ Create relevant plots for the hvarma order finder algorithm. """
    ps = []
    pos_freq = []
    neg_freq = []
    for p in orders:
        ps.append(p)
        freqs = orders[p].get_frequency(20)
        pos_freq.append(freqs[0])
        neg_freq.append(freqs[2])

    ps = np.array(ps)
    pos_freq = np.array(pos_freq)
    neg_freq = np.array(neg_freq)

    fig1 = pretty_plot_search(ps, pos_freq, neg_freq, found_p, stat_name)
    fig2 = pretty_plot_search_errors(ps, pos_freq, neg_freq, found_p, stat_name, tol)

    if output_dir is not None:
        filename1 = os.path.join(output_dir, f'{stat_name}_order_search_freq.png')
        filename2 = os.path.join(output_dir, f'{stat_name}_order_search_diff.png')

        fig1.savefig(filename1, dpi=300, bbox_inches='tight')
        fig2.savefig(filename2, dpi=300, bbox_inches='tight')

    return fig1, fig2


def pretty_plot_search(ps, pos_freq, neg_freq, found_p, stat_name, y_label='Frequency', tol=None):
    """ Plot the successive frequency values and highlight an order of choice """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    # Set Times font
    plt.rc('font', family='serif', size=12)
    plt.rcParams['mathtext.fontset'] = 'dejavuserif'

    fig = plt.figure()
    plt.rc('axes', axisbelow=True)
    plt.grid()
    plt.plot([min(ps), max(ps)], [0, 0], 'k-', linewidth=1)

    plt.scatter(ps, pos_freq, label='Positive', facecolors='none', edgecolors='k')
    plt.scatter(ps, neg_freq, label='Negative', facecolors='k', edgecolors='k')
    plt.scatter(found_p, pos_freq[np.argmax(ps == found_p)], facecolors='none', edgecolors='r')
    plt.scatter(found_p, neg_freq[np.argmax(ps == found_p)], facecolors='r', edgecolors='r')

    if tol is not None:
        plt.plot([min(ps), max(ps)], [tol, tol], 'k--', linewidth=1)
        plt.yscale('log')

    plt.minorticks_on()
    plt.xlabel('Model order')
    plt.ylabel(y_label)
    plt.text(0.85 * max(ps) + 0.15 * min(ps), 0.88*max(max(pos_freq), max(neg_freq)),
             stat_name, fontweight='bold', fontsize=14)
    return fig


def pretty_plot_search_errors(ps, pos_freq, neg_freq, found_p, stat_name, tol=0.1):
    """ Plot the absolute value of the frequency differences in consecutive model orders. """
    ids = ps.argsort()
    ps = ps[ids]
    pos_freq = pos_freq[ids]
    neg_freq = neg_freq[ids]
    diffs_pos = []
    diffs_neg = []
    orders = []
    for p in ps:
        if p-3 in ps:
            orders.append(p)
            pos = abs(pos_freq[np.argmax(ps == p)] - pos_freq[ps == p-3])
            neg = abs(neg_freq[np.argmax(ps == p)] - neg_freq[ps == p-3])
            diffs_pos.append(neg)
            diffs_neg.append(pos)

    return pretty_plot_search(np.array(orders), np.array(diffs_pos), np.array(diffs_neg),
                              found_p, stat_name, 'Freq. diff', tol=tol)
