# Rodionov_STARS
A command line implementation of the Student's t-test for the Analysis of Regime Shifts (STARS)

The STARS algorithm aims to identify regime changes in time series data. It is described in detail by Rodionov (2004). Many implementations of this algorithm are available online, notably from the author in the form of a [Microsoft Excel plugin](https://www.beringclimate.noaa.gov/regimes/). The objective of this one, is to provide a tool that can be used by collaborators familiar with different programming langagues, as long as they are familiar with the Linux shell and have a basic Python installation. For Python users, it can also be imported into scripts as a module.

## Installation

For users without a Python installation, a convenient way to install all the needed dependencies is install via [Anaconda](https://www.anaconda.com/download). Otherwise, ensure that the following modules are installed to your system (or virtual environment):

* Scipy
* Numpy
* Pandas
* Matplotlib
* Argparse
* os

Once Python is installed, clone this repository into your Python path. Typically, with an Anaconda install, this would be done in one of the following locations:

```python
/path/to/Anaconda/lib/pythonX.X/site-packages/
/path/to/Anaconda/envs/environment_name/lib/pythonX.X/site-packages/
```

To test the installation, you can run the validation script, which reproduces the results from the algorithm example shown in Rodionov (2004) using the Pacific Decadal Oscillation (PDO) January index.

```python
$ python rodionov_stars_validation.py
```



## References

Rodionov, S.N. (2004) A sequential algorithm for testing climate regime shifts. Geophys. Res. Lett., 31, L09204, doi:10.1029/2004GL019448.
