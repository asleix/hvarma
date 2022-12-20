import sys
import time
import matplotlib.pyplot as plt
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
    station_name = average_data.station
    num_windows = average_data.num_windows
    if write_png:
        if param.oname == 'default':
            outfile = 'output/{}_p{}_win{}'.format(station_name, param.p, num_windows)
        else:
            if not param.oname.endswith('/'):
                param.oname += '/'
            outfile = param.oname + '{}_p{}_win{}'.format(station_name, param.p, num_windows)
        outfile += '.png'
    else:
        outfile = None

    err = param.plot_conf/2
    pos_freq, pos_err, neg_freq, neg_err = average_data.get_frequency(param.freq_conf)

    return pretty_show(param.f0, param.f1, param.npun, average_data.get_spectrum_percentile(50),
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
