[![Python package](https://github.com/asleix/hvarma/actions/workflows/python-package.yml/badge.svg)](https://github.com/asleix/hvarma/actions/workflows/python-package.yml)

# H/V-ARMA calculator

The horizontal to vertical (H/V) spectral ratio is used to estimate
the resonant frequency of sediments. HV-ARMA does so by 
estimating the transfer function in a surface layer for
three-dimensional micro-tremor seismogram data.

This HV-ARMA module and the provided example scripts allow
to find the resonant frequency in micro-tremor seismograms
for a given HV-ARMA model order. 

- `hvarma`: Autoregressive moving average model of microtremor H/V spectral ratio. Software Impacts. 2025. [[Link to Paper](https://doi.org/10.1016/j.simpa.2025.100745)].

![Resonant frequencies](https://raw.githubusercontent.com/asleix/hvarma/main/examples/b001_p42.png)

In addition, we provide 
guidance for finding the best model order for given input 
data.

## Features

This software allows to compute:
- Spectral ratio, with lower and upper error bounds.
- Estimated resonant frequency.
- Signal coherence.
- Optimal model order.

The supported data format by default are [SAC](https://ds.iris.edu/ds/support/faq/17/sac-file-format/) binaries (Seismic Analysis 
Code). However, the 
module API accepts any float-valued numpy array.

## Installation guide

This software combines Python and C code to speed up 
critical sections. Hence, Python 3 and `gcc` installations 
are required.

Using a virtual environment is recommended if 
you have multiple projects. Instructions are added below. If 
you choose not to use it, the following instructions
might need to replace python and pip with `python3` and `pip3`.

Install the python module. 
From the root of the H/V ARMA folder, use

```
pip install .
```

Now you can use `hvarma` by simply importing the module in your
python scripts using `import hvarma`. You can also use the
default scripts such as `run.py`. Feel free to move `run.py`
to your working directory.

### Using a virtual environment

You should first install `virtualenv` using
``` 
pip3 install virtualenv
```
Then, to create a virtual environment, use
``` 
python3 -m virtualenv venv
```
This will create the folder `venv` in your current directory.
You can place it in your working directory. To activate
the virtual environment (you will need to do it once per 
session)
``` 
source venv/bin/activate
```
Now that your environment is active, notice that your
python3 is now `python` and your pip is `pip`, as you
can verify with
```
which python; which pip
```


## Using the example scripts

To calculate the horizontal-to-vertical ratio on your data
using default parameters, use

```
python run.py Z_data_filename.sac N_data_filename.sac E_data_filename.sac
```

You can specify the model parameters in a file. In the examples,
the file `args.txt` specifies different parameters. You can
indicate the file with the argument `--args=args.txt` when
running `python run.py`. To get a list with all the 
possible arguments, please run
```
python run.py -h
```

The output is stored by default in the current directory unless
specified with the argument `--output_dir=output/`. If an argument
is specified both in `args` and as a command line argument,
the latter will override the other values.

We implement an algorithm to find a candidate for a 
good model order to describe the signal. The script is 
`find_model_order.py` and can be used in the same
fashion as `run.py`.


## Using the module

The module `hvarma` implements different functions and
data structures.

- `Data`. Class to read SAC file data.
- `ArmaParam`. Class to store the model parameters.
- `run_model`. Run a model specified by an ArmaParam instance on an
                instance of Data.
- `plot_hvarma`. Plot your results from run_model in a plot
                like the ones show in this readme.
- `find_optimal_order`. Execute a fast algorithm to 
              test different candidate model orders 
              and choose the smallest that satisfies 
              a convergence condition.


## Parameter specification

Parameters are defined in the args.txt file.
- `model_order`: Number of lags "p" of the ARMA(p, p) model.
- `maxtau`: Number of considered lags in the covariances. The value should be at most `window_size/2`.
- `mu`: Weight of the forward prediction error in coefficient estimation.
- `nu`: Weight of the backward prediction error in coefficient estimation.
- `nfir`: Number of covariance lags considered in the correlation 
  matrices during the coherence calculation. A relatively small number
  should suffice.
- `neg_freq`: Starting point of the frequency interval.
- `pos_freq`: Ending point of the frequency interval
- `freq_points`: Number of points within the frequency interval at which
  to compute the spectral ratio.
- `window_size`: Number of data points for each time window.
- `overlap`: Number of overlapping points between two consecutive windows.
- `max_windows`: Maximum number of windows. Data beyond the max window will be ignored.
  If max_windows is higher than available data, max_windows will be ignored.
- `freq_conf`: Confidence interval for the resonant frequency. 
  The lower the number, higher the confidence. It is in the interval 0-100. 
- `plot_conf`: Confidence interval of the spectral ratio.
- `output_dir`: Names will be [stationname]_p[arma_order]_win[number_of_windows].
               If output_path is 'default', files will be stored in .../output 


The default arguments are (as in an `args.txt`)

```
model_order=74
maxtau=128
mu=0.5
nu=0.5
nfir=40
neg_freq=-20
pos_freq=20
freq_points=1024
window_size=512
overlap=256
max_windows=1000
freq_conf=20
plot_conf=50
output_dir=.
```






