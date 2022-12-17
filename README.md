# H/V-ARMA calculator

The H/V method is the horizontal to vertical spectral ratio 
and it is used to estimate the resonant frequency of sediments.
This H/V-ARMA calculator allows for computing this spectral ratio 
from seismographic data saved in SAC format.

## Features

This software allows to compute:
- Spectral ratio, with lower and upper error bounds.
- Estimated resonant frequency.
- Signal coherence.

Input: 3 SAC data files and a parameters file, "args.txt".
Output: Text file with spectral ratio and coherence, and a plot.

## Installation guide

This software combines Python and C code to speed up 
critical sections. Thus, both Python 3 and C languages
must be installed on the computer. This software will not work
in Python 2. 

Install python dependencies. From the root of the H/V ARMA folder, use

```
pip3 install -r requirements.txt
```

Compile C shared library via Makefile. Compile using

```
make
```

After these steps, the file obj/gradient.so should be created and
everything is set up.

In order to recompile C library, use 
``` make clean ``` first.


## Usage

After setting up args.txt file, run

```
python3 run.py Z_data_filename.sac N_data_filename.sac E_data_filename.sac
```

Output is stored by default in either /output/ or the specified *output_name* parameter in args.txt.


## Parameter specification

Parameters are defined in the args.txt file.
- arma_order: Number of lags "p" of the ARMA(p, p) model.
- maxtau: Number of considered lags in the covariances. At most,
  it should be window_size/2.
- forward_weight: Weight of the forward prediction error in coefficient estimation. Set to 0 or "sigma" to use vertical variance.
- backward_weight: Weight of the forward prediction error in coefficient estimation. Set to 0 or "sigma" to use horizontal variance.
- nfir: Number of covariance lags considered in the correlation 
  matrices during the coherence calculation. A relatively small number
  should suffice.
- ini_freq: Starting point of the frequency interval.
- fin_freq: Ending point of the frequency interval
- num_points: Number of points within the frequency interval at which
  to compute the spectral ratio.
- window_size: Number of data points for each time window.
- window_shift: Number of overlapping points between two consecutive windows.
- max_windows: Maximum number of windows. Data beyond the max window will be ignored. 
  If max_windows is higher than available data, max_windows will be ignored.
- freq_conf_interval: Confidence interval for the resonant frequency. 
  The lower the number, higher the confidence. It is in the interval 0-100. 
- spectral_conf_interval: Confidence interval of the spectral ratio.
- output_path: Names will be [stationname]_p[arma_order]_win[number_of_windows].
               If output_path is 'default', files will be stored in .../output 


The file args.txt should look like

```
arma_order=74
maxtau=256
forward_weight=0.5
backward_weight=0.5
nfir=40
ini_freq=-20
fin_freq=20
num_points=512
window_size=512
window_shift=256
max_windows=1000
freq_conf_interval=20
spectral_conf_interval=50
output_name=default
```

## Program structure

- /run.py. Main script that handles workflow. It is the one to be called.
- /args.txt. Arguments file. Filename should not be changed.
- /include/read_input.py. Contains data and parameter reading functions.
- /include/write_output.py. Display and saving functions, including 
                           a progressbar, plot and ascii text saving.
- /include/processing.py. Handle window processing and all-window averaging.
- /include/compute.py. Hard-working formula functions.
- /src/gradient.c. Most expensive computation, the optimality conditions matrix. 
- /bin/gradient.so. Shared library, created at C compile time.
- /output/. Default output directory. Output is stored here if parameter output_path=default.

We also provide additional data files and output as an example computation.
They can be found in data and output folders.

File tree:
```
.
├── Makefile
├── README.md
├── args.txt
├── data
│   ├── B001_E.sac
│   ├── B001_N.sac
│   └── B001_Z.sac
├── include
│   ├── compute.py
│   ├── processing.py
│   ├── read_input.py
│   └── write_output.py
├── obj
│   └── gradient.so
├── output
│   ├── out.png
│   └── out.txt
├── requirements.txt
├── run.py
├── src
│   └── gradient.c
```




