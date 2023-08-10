#!/Users/roy/anaconda3/bin/python

"""

This script is used to calculate filament age and filament depolymerization velocity.
Filament age is calculated as a function of time spent between polymerization and depolymerization.
Depolymerization velocity is calculated with varying window sizes.
Fit parameters for a kinetic model of 1/v_depol as a function of filament age are calculated for different window sizes.

Exponential decay model is expressed as:
y = a*exp(-bx) + c
where:
x = filament age
y = 1/v_depol
a + c = 1/v_ADP-Pi_depol
b = kr
c = 1/v_ADP_depol

Input: <cumulative-filament-length-data-file.csv>
Output:
    <filename>_velocity_all-datasets.csv
    <filename>_velocity_optimal-datasets.csv
    <filename>_velocity_all-parameters.csv
    <filename>_velocity_optimal-parameters.csv

Input data file should contain the following columns:
norm_time       -> time since depolymerization (s)
mean_monomers   -> mean filament length in monomers calculated from multiple single filament tracked data
sd [optional]   -> standard deviation of mean filament length
se [optional]   -> standard error of mean filament length

Parameters for the exponential decay model are estimated for each window size.
Optimal parameters are obtained such that the average percentage standard deviation error in estimation is less than a threshold value (10% recommended).

Velocity dataset files contain the following fields:
norm_time       -> as in input file
mean_monomers   -> as in input file
sd              -> as in input file
se              -> as in input file
age             -> filament age (s)
velocity        -> polymerization velocity (monomers/s)
norm_velocity   -> depolymerization velocity 1/v_depol (s/monomer)
window_size     -> window size used for calculating velocity of the centered time point

Parameter files contain the following fields:
window_size     -> window size used for calculating velocity of the centered time point
a               -> predicted fit parameter a
b               -> predicted fit parameter b
c               -> predicted fit parameter c
a_err           -> standard deviation error in prediction of a
b_err           -> standard deviation error in prediction of b
c_err           -> standard deviation error in prediction of c
a_per_err       -> percentage standard deviation error in prediction of a
b_per_err       -> percentage standard deviation error in prediction of b
c_per_err       -> percentage standard deviation error in prediction of c
avg_per_err     -> average percentage standard deviation error in parameter prediction
sym_window_size -> symmetric window unit size

"""

__author__ = "Ankit Roy"
__copyright__ = "Copyright 2023, Bieling Lab, Max Planck Institute of Molecular Physiology"
__credits__ = ["Ankit Roy"]
__license__ = "GPLv3"
__maintainer__ = "Ankit Roy"
__email__ = "ankitroy.post@gmail.com"
__status__ = "Production"
__version__ = "1.0"

import numpy as np
import pandas as pd
import sys
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import seaborn as sns

# Get data
def dataIN(filename):
    data = pd.read_csv(filename)
    return data

# Subset polymerization data
def get_pol_data(data, min_norm_time):
    return data.loc[(data.norm_time <= 0) & (data.norm_time >= min_norm_time),].copy()

# Subset depolymerization data
def get_depol_data(data, max_norm_time):
    return data.loc[(data.norm_time >= 0) & (data.norm_time <= max_norm_time),].copy()

# Linear fit function
def local_LinFit(x, y):
    m, b = np.polyfit(x, y, 1)
    return m

# Calculate filament age
def calc_age(data, pol_velocity):
    return np.abs(data.mean_monomers/pol_velocity) + data.norm_time

# Calculate instant depolymerization velocity
def calc_instantVelocity(time, monomers, size=5):
    velocities = {}                 # stores velocities calculated over time windows

    # calculate instant velocity from a sliding window
    for n in range(len(time) - size + 1):
        x = time[n:n+size]          # time window
        y = monomers[n:n+size]      # monomer window

        v = local_LinFit(x, y)      # velocity over the specified window
        t = time[n+int(size/2)]     # time point around which the window was applied

        velocities[t] = v           # store time and velocity pair

    return velocities

# Calculate inverse of depolymerization velocity
def calc_inverse_velocity(data):
    return 1/np.abs(data.velocity)

# Exponential decay fit function
def expDecay(x, a, b, c):
    return a * np.exp(-x * b) + c

# Clean data
def clean_data(data, max_norm_time):
    data_clean = data.dropna()
    data_clean =  data_clean.loc[data_clean.norm_time <= max_norm_time,]
    return data_clean

# Calculate optimum fit parameters
def get_fit_parameters(fit_function, data, a_bound=1, b_bound=0.1, c_bound=0.21):
    popt, pcov = curve_fit(fit_function,
                           data.age,
                           data.norm_velocity,
                           bounds=(0, [a_bound, b_bound, c_bound])
                           )
    
    perr = np.sqrt(np.diag(pcov))

    return popt, perr

# Fit data with a specific window size
def fit_data_single(filename, depol_data, max_norm_time, window_size, a_bound, b_bound, c_bound):

    # Calculate instantaneous depolymerization velocity
    depol_velocity = calc_instantVelocity(list(depol_data.age),
                                          list(depol_data.mean_monomers),
                                          window_size)
    depol_data['velocity'] = depol_data['age'].map(depol_velocity)

    # Calculate the inverse of depolymerization velocity
    depol_data["norm_velocity"] = calc_inverse_velocity(depol_data)

    # Clean up
    depol_data_clean = clean_data(depol_data, max_norm_time)

    # Add a window size identifier
    depol_data_clean["window_size"] = window_size

    # Optimize fit parameters
    popt, perr = get_fit_parameters(expDecay, depol_data_clean, a_bound, b_bound, c_bound)

    return popt, perr, depol_data_clean

# Fit data with multiple window sizes
def fit_data_multiple(filename, depol_data, max_norm_time, window_sizes, a_bound, b_bound, c_bound):

    pstats = {}             # stores optimum parameters and their errors for all window sizes
    depol_datasets = {}     # store depolymerization data for all datasets with velocity calculated using various window sizes

    # Find optimum fit parameters for all window sizes
    for win_size in window_sizes:
        popt, perr, depol_data_clean = fit_data_single(filename,
                                                       depol_data,
                                                       max_norm_time,
                                                       win_size,
                                                       a_bound,
                                                       b_bound,
                                                       c_bound)
        
        # Store parameter stats
        pstats[win_size] = np.concatenate((popt, perr), axis=0)

        # Store depolymerization data
        depol_datasets[win_size] = depol_data_clean

    # Construct parameter statistics dataframe
    pstats_df = pd.DataFrame.from_dict(pstats,
                                  orient='index',
                                  columns=['a', 'b', 'c', 'a_err', 'b_err', 'c_err'])
    pstats_df.index.name = "window_size"
    
    pstats_df["a_per_err"], pstats_df["b_per_err"], pstats_df["c_per_err"] = [
        pstats_df.a_err/pstats_df.a*100,
        pstats_df.b_err/pstats_df.b*100,
        pstats_df.c_err/pstats_df.c*100]
    
    pstats_df["avg_per_err"] = (pstats_df.a_per_err + pstats_df.b_per_err + pstats_df.c_per_err)/3
    
    return pstats_df, depol_datasets

# Get optimal parameter fits
def get_optimal_fits(parameters, avg_per_err_threshold=20):
    return parameters.loc[parameters.avg_per_err <= avg_per_err_threshold]

# Combine datasets analyzed at different window sizes
def combine_datasets(depol_datasets):
    datasets_combined = pd.concat(list(depol_datasets.values()), ignore_index=True)
    return datasets_combined

# Get datasets where velocities were calculated with optimal parameters
def get_optimal_datasets(datasets_combined, parameters_optimal):
    
    # Optimal window sizes
    optimal_window_sizes = list(parameters_optimal.index)
    
    # Subset of datasets that were analyzed with optimal window sizes
    optimal_datasets = datasets_combined.loc[
        datasets_combined.window_size.isin(optimal_window_sizes), 
        ]
    
    return optimal_datasets

# Calculate symmetrical unit window size
def get_window_sym_unit(data):
    data_copy = data.copy()
    data_copy["sym_window_size"] = (data_copy.index - 1)/2
    return data_copy

# Plot data points and fit function
def plot_data(datasets, fit_function, popt):

    sns.lineplot(
        data = datasets,
        x = datasets.age,
        y = datasets.norm_velocity,
        marker = 'o',
        dashes = False,
        errorbar = 'sd',
        label = 'data',
        color = 'black'
    )

    sns.lineplot(
        data = datasets,
        x = datasets.age,
        y = fit_function(datasets.age, *popt),
        label = 'fit',
        color = 'red'
    )

    plt.xlabel(r"Filament age $\tau$ [s]")
    plt.ylabel(r"$1/v_{depol}$ [s/monomers]")
    plt.title("Instantaneous Depolymerization Velocity")
    # plt.xlim(0, 430)
    # plt.ylim(0.05,0.4)

    plt.show()

# Plot error in parameter prediction for different window sizes
def plot_error(toggle, parameters, parameters_optimal, opt_param_threshold):

    min_opt_win_size = parameters_optimal.sym_window_size.min()
    max_opt_win_size = parameters_optimal.sym_window_size.max()

    sns.lineplot(
        data = parameters,
        x = parameters.sym_window_size,
        y = parameters.avg_per_err,
        marker = 'o',
    )

    plt.axhline(opt_param_threshold, ls='--', c='black')
    plt.axvline(min_opt_win_size, ls='--', c='black')
    plt.axvline(max_opt_win_size, ls='--', c='black')
    plt.ylim(0, 50)
    plt.xlabel("Window Symmetric Unit Size")
    plt.ylabel("Average Percentage Standard Deviation Error")
    plt.title("Fit Parameter Prediction Error")

    plt.show()

# Write output file
def dataOUT(filename, data, suffix, index=False):
    outname = '.'.join(filename.split('.')[:-1])
    outname = f'{outname}_velocity_{suffix}.csv'
    data.to_csv(outname, index = index)

# Write ouput files when asked
def write_output_files(toggle, filename, datasets_combined, optimal_datasets, parameters, parameters_optimal):

    # Write output files if toggled to True
    if toggle:
        dataOUT(filename, datasets_combined, "all-datasets")
        dataOUT(filename, optimal_datasets, "optimal-datasets")
        dataOUT(filename, parameters, "all-parameters", True)
        dataOUT(filename, parameters_optimal, "optimal-parameters", True)
        print("Output files written.")
    else:
        print("Output files suppressed.")

# Main function
def main():
    filename = sys.argv[1]                  # input file name
    window_sizes = np.arange(3, 71, 2)      # window sizes to iterate over
    min_norm_time = -150                    # minimum time from depolymerization start
    max_norm_time = 400                     # maximum time from depolymerization start
    a_bound = 1                             # upper bound for estimation of a
    b_bound = 0.1                           # upper bound for estimation of b
    c_bound = 0.21                          # upper bound for estimation of c
    opt_param_threshold = 10                # maximum average percentage standard deviation error in fit parameters
    toggle_output_writing = True            # Set True to write output files
    toggle_error_plot = True                # Set True to show error plot

    # Get input data
    data = dataIN(filename)

    # Subset polymerization and depolymerization data
    pol_data = get_pol_data(data, min_norm_time)
    depol_data = get_depol_data(data, max_norm_time)

    # Calculate polymerization velocity
    pol_velocity = local_LinFit(pol_data.norm_time, pol_data.mean_monomers)

    # Calculate filament age
    depol_data["age"] = calc_age(depol_data, pol_velocity)

    # Fit data
    parameters, depol_datasets = fit_data_multiple(filename,
                                                   depol_data,
                                                   max_norm_time,
                                                   window_sizes,
                                                   a_bound,
                                                   b_bound,
                                                   c_bound)
    
    # Combine datasets
    datasets_combined = combine_datasets(depol_datasets)

    # Get a subset of optimal parameters
    parameters_optimal = get_optimal_fits(parameters, opt_param_threshold)

    # Get symmetrical unit size for windows
    parameters = get_window_sym_unit(parameters)
    parameters_optimal = get_window_sym_unit(parameters_optimal)

    # Get optimal datasets
    optimal_datasets = get_optimal_datasets(datasets_combined, parameters_optimal)

    # Mean of optimal parameters
    a_opt_mean = parameters_optimal.a.mean()
    b_opt_mean = parameters_optimal.b.mean()
    c_opt_mean = parameters_optimal.c.mean()

    # Show optimal parameters
    print(f"Average of optimal parameters [a, b, c]: [{a_opt_mean:.4f}, {b_opt_mean:.4f}, {c_opt_mean:.4f}]")

    # Write output files
    write_output_files(toggle_output_writing, filename, datasets_combined, optimal_datasets, parameters, parameters_optimal)

    # Plot data with fit
    plot_data(optimal_datasets, expDecay, [a_opt_mean, b_opt_mean, c_opt_mean])

    # Plot parameter error over different window sizes
    plot_error(toggle_error_plot, parameters, parameters_optimal, opt_param_threshold)

main()

# Ankit Roy
# 25th May, 2023
#   >> Updated script to calculate filament age and fit velocity vs filament age.
# 26th May, 2023
#   >> Now added feature to calculate fit parameters with instantaneous velocity estimations for multiple window sizes.
#   >> Generates a combined dataset with velocities calculated from different window sizes.
#   >> Generates a file with fit parameters and its standard deviation errors.
#   >> Allows filtering optimal parameters based on the average standard deviation error of fit parameters.
#   >> Generates subset files with optimal parameters and velocities calculated with optimal window sizes.
#   >> Plots an aggreagate of 1/v_depol vs filament age with a fit obtained from the mean of optimal parameters.
#   >> Generates 4 output files:
#           - <file-name>_all-datasets.csv : velocity data calculated for all window sizes
#           - <file-name>_optimal-datasets.csv : velocity data calculated for optimal window sizes
#           - <file-name>_all-parameters.csv : parameter stats calculated for all window sizes
#           - <file-name>_optimal-parameters.csv : parameter stats calculated for optimal window sizes
#   >> Added to toggle to chose if output files are wanted to be written.
#   >> Generates a plot of parameter prediction average error when toggle is switched to True.
# 10th August, 2023
#   >> Updated docstrings.
