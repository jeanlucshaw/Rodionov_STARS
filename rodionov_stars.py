"""
Command-line implementation of the Sequential t-test Analysis of Regime Shifts (STARS)

User instructions:

    This script can be used from command line in any linux shells. Ensure that your
    Python environment has the Pandas, Scipy, Numpy, Matplotlib and Argparse modules
    installed. A detailed description of the option flags in this use case can be
    accessed by running:

        $ python rodionov_stars.py -h

    To illustrate command line usage, regimes from the example in the original STARS
    paper (Rodionov, 2004) can be identified by running:

        $ python rodionov_stars.py data/rodionov_2004_figure_2_PDO.csv example_output.csv -t 0 -c 1

    This script can also be used as a module and its functions imported into new Python
    scripts, assuming that the containing folder and its contents have been placed on
    the user's Python path.

Input data format:

    The input data should be supplied as a comma-separated ASCII column file (.csv)
    with no header. See `data/rodionov_2004_figure_2_PDO.csv` for an example.

References:

    Rodionov, Sergei. (2004) A sequential algorithm for testing climate regime shifts
    Geophysical Research Letters. 31(9). 1--5

Differences with the algorithm described in the paper and precisions:

    1- The search for regime shift starts at index 1, instead of `l + 1` (0, 1, 2, ...).
    2- The regime shift mean is not updated unless the index is more than `l` points away
       from the last change point.
    3- the variance used to calculate `diff` is the biased/sample/divide-by-N variance.

Assumptions about the data:

    1- data are spaced at regular time intervals
    2- no missing data (bad values or gaps)

For questions and bug reports, contact Jean-Luc.Shaw@dfo-mpo.gc.ca

"""
from scipy.stats import t as t_distribution
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse as ap
import os


def calculate_rsi(variable, index, x_star, l, sigma_l):
    """
    Regime change index (RSI) calculations

    Parameters
    ----------
    variable : pandas.Series
        the data used to detect regime shifts
    index : int
        of the change point value for which to calculate the RSI
    x_star : func
        how to calculate x_star; either (x_i - x_R2) or (x_R2 - x_i)
    l : int
        cut-off length of the regimes to be determined
    sigma_l : float
        average of the rolling (l-length) variance of `variable`.

    Returns
    -------
    rsi : float
        the RSI value for `variable` at `index`

    """
    rsi = 0
    for j_ in np.arange(index, min(index + l, variable.size)):
        rsi += x_star(variable.iloc[j_]) / l / np.sqrt(sigma_l)

        # Step 6: no regime change
        if rsi < 0:
            rsi = 0
            break

    return rsi


def plot(filename, rodionov_results):
    """
    Make a quick plot of the STARS results

    Mostly intended for use via the command line (option -P) to get a quick
    idea of the analysis results. For more tailored results, users can work
    from the `csv` output data file.

    Parameters
    ----------
    filename : str
        under which to save the results figure.
    rodionov_results : pandas.DataFrame
        containing the RSI and regime_id output by `regime_shift_rodionov`, and the
        time series data `variable`. The time axis `time` will be used if supplied.

    """
    # Init
    _, ax = plt.subplots(2, figsize=(10, 6), sharex=True)

    # Manage time input presence/absence
    if 'time' in rodionov_results.keys():
        time = rodionov_results.time
    else:
        time = np.arange(rodionov_results.RSI.size)

    # ========
    # Variable data

    # draw the detected regimes as alternating colored bars
    rbcm = np.array(['lightskyblue', 'lightcoral'])
    colors = rbcm[np.int64(rodionov_results.regime_id.values % 2)]
    ax[0].bar(time, rodionov_results.variable, color=colors)

    # ========
    # RSI data
    index_change = rodionov_results.RSI > 0
    ax[1].bar(time[index_change], rodionov_results.RSI[index_change], color='gray')

    # Plot parameters
    for a_ in ax:
        a_.tick_params(which='both', top=False, right=False)
        a_.tick_params(which='minor', bottom=False, left=False)

    # Save plot next to output data
    plt.savefig(filename, dpi=300)
    
    return None


def regime_shift_rodionov(variable, l=10, p=0.05):
    """
    Apply the Rodionov (2004) STARS algorithm to detect regime shifts

    Parameters
    ----------
    variable : 1D numpy.ndarray or pandas.Series
        the time series in which to detect regime shifts
    l : int
        the cut-off length of the regimes to be determined
    p : float
        probability used in the Student's t-test.

    Returns
    -------
    rsi : 1D numpy.ndarray of float
        the regime shift index calculated at each time step
    regime_id : 1D numpy.ndarray
        labels attached to each time step assigning them to a regime
    regime_mean : 1D numpy.ndarray of float
        used by the algorithm at each time step
    diff : float
        parameter which determines the regime change thresholds around `regime_mean`

    """

    # Ensure input is a Pandas series
    variable = pd.Series(variable)

    # Initialize arrays
    rsi = np.zeros(variable.size)
    regime_id = np.zeros(variable.size)
    regime_mean = np.zeros(variable.size) * np.nan
    regime_count = 0
    change_point = 0

    # Step 2: calculate the diff parameter
    sigma_l = (variable.rolling(window=l).var(ddof=0)).mean()  # note: biased/sample variance
    diff = t_distribution.ppf(1 - p/2, 2*l - 2) * np.sqrt(2 * sigma_l / l)

    # Step 3: calculate the initial regime and regime change thresholds
    x_bar_R1 = variable.iloc[:l].mean()

    # Step 4: iterate over new values and test for a regime change
    for i_ in np.arange(1, variable.size):

        # Save current regime mean value
        regime_mean[i_] = x_bar_R1

        # Possible regime change?
        if variable.iloc[i_] > (x_bar_R1 + diff):

            # Step 5
            rsi[i_] = calculate_rsi(variable,
                                    i_,
                                    lambda x: x - (x_bar_R1 + diff),
                                    l,
                                    sigma_l)

        elif variable.iloc[i_] < (x_bar_R1 - diff):

            # Step 5
            rsi[i_] = calculate_rsi(variable,
                                    i_,
                                    lambda x: (x_bar_R1 - diff) - x,
                                    l,
                                    sigma_l)
        else:
            rsi[i_] = 0

        # Step 7: Regime change
        if rsi[i_] > 0:
            x_bar_R1 = variable.iloc[i_ : min(i_ + l, variable.size)].mean()
            regime_count += 1
            change_point = i_

        # Step 6: No regime change
        else:
            if (i_ - change_point) > l:
                x_bar_R1 = variable.iloc[i_ - l + 1: i_ + 1].mean()

        # Assign current regime ID to this value
        regime_id[i_] = regime_count

    return rsi, regime_id, regime_mean, diff


# ======================
# Command line interface
# ======================
if __name__ == '__main__':
    parser = ap.ArgumentParser(prog=__doc__)
    parser.add_argument('inputfile',
                        metavar='',
                        help='csv file containing the time series')
    parser.add_argument('outputfile',
                        metavar='',
                        help='csv file containing STARS results (RSI, REGIME ID, REGIME MEAN, DATA, [TIME])')
    parser.add_argument('-l', '--cutoff-length',
                        metavar='',
                        type=int,
                        default=10,
                        help='cutoff length of the regimes to be determined (default: 10)')
    parser.add_argument('-p', '--probability',
                        metavar='',
                        type=float,
                        default=0.05,
                        help="probability used in the Student's t-test (default: 0.05)")
    parser.add_argument('-c', '--data-column',
                        metavar='',
                        type=int,
                        default=0,
                        help='the column of the csv file to use as data (defaults: leftmost)')
    parser.add_argument('-t', '--time-column',
                        metavar='',
                        type=int,
                        default=None,
                        help='the column of the csv file containing the time axis (defaults: not used)')
    parser.add_argument('-P', '--plot',
                        # metavar='',
                        action='store_true',
                        help='Save a plot of the identified regimes and change point RSIs')
    args = parser.parse_args()

    # Read input data
    variable = pd.read_csv(args.inputfile, header=None).iloc[:, args.data_column]
    if args.time_column != None:
        time = pd.read_csv(args.inputfile, header=None).iloc[:, args.time_column]
        
    # Run the algorithm
    rsi, regime_id, regime_mean, diff = regime_shift_rodionov(variable,
                                                              l=args.cutoff_length,
                                                              p=args.probability
                                                              )

    # Format and output
    output = pd.DataFrame({'RSI': rsi,
                           'regime_id': regime_id,
                           'regime_mean': regime_mean,
                           'variable': variable})
    if args.time_column != None:
        output.loc[:, 'time'] = time
    output.to_csv(args.outputfile, index=False)

    # Optional plot
    if args.plot:
        figure_name, _ = os.path.splitext(args.outputfile)
        plot(f'{figure_name}.png', output)
