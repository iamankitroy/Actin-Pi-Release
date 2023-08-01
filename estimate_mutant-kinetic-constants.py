#!/Users/roy/anaconda3/bin/python

__author__ = "Ankit Roy"
__copyright__ = "Copyright 2023, Bieling Lab, Max Planck Institute of Molecular Physiology"
__credits__ = ["Ankit Roy"]
__license__ = "GPLv3"
__maintainer__ = "Ankit Roy"
__email__ = "ankitroy.post@gmail.com"
__status__ = "Production"
__version__ = "1.0"

import sys
import pandas as pd
import numpy as np

# Get data
def dataIN(filename):
    return pd.read_csv(filename)

# Calculate average for parameter c of mutant
def get_average_c(data):
    return data.norm_velocity.mean()

# Calculate parameter a for mutant
def get_a_mut(a_wt, c_wt, c_mut, multiplier=1):
    # 1/v_adp_pi_depol = a_wt + c_wt
    # mutant 1/v_adp_pi_depol can be scaled differently to wt using the multiplier
    a_mut = (a_wt + c_wt)*multiplier - c_mut
    return a_mut

# Exponential decay function with parameters a, b and c
def expDecay(x, a, b, c):
    return a * np.exp(-b * x) + c

# Cost function
def cost_function(y_pred, y):
    return y_pred - y

# Optimisation function
def optimize(b_min, b_max, x, y, a, c, tolerance):
    b_all = np.arange(b_min, b_max, 0.001)      # parameter b search space
    b_optimum = b_min                           # optimum b

    # search for optimum b
    for b_i in b_all:
        y_pred = expDecay(x, a, b_i, c)         # y prediction using parameter b_i
        cost = cost_function(y_pred, y)         # cost of prediction
        
        # update optimum b when cost is tolerated
        if cost <= tolerance:
            b_optimum = b_i
            break

    return b_optimum

# Main function
def main():
    filename = sys.argv[1]                      # mutant data file with normalised velocity
    data = dataIN(filename)                     # mutant data

    a_wt = 0.18350639                           # wt exponential decay parameter a
    c_wt = 0.15155544                           # wt exponential decay parameter c
    tolerance = 1e-3                            # cost tolerance
    
    c_mut = get_average_c(data)                 # mutant exponential decay parameter c

    # mutant exponential decay parameter a; assuming mutant v_adp_pi_depol = wt
    a_mut_1 = get_a_mut(a_wt, c_wt, c_mut, 1)
    # mutant exponential decay parameter a; assuming mutant v_adp_pi_depol = 2x wt
    a_mut_2 = get_a_mut(a_wt, c_wt, c_mut, 0.5)

    x = data.norm_time[0]                       # first instance of norm_time
    y = data.norm_velocity[0]                   # corresponding norm_velocity to x

    # mutant exponential decay parameter b; assuming mutant v_adp_pi_depol = wt
    b_mut_1 = optimize(0.010, 0.300, x, y, a_mut_1, c_mut, tolerance)
    # mutant exponential decay parameter b; assuming mutant v_adp_pi_depol = 2x wt
    b_mut_2 = optimize(0.010, 0.300, x, y, a_mut_2, c_mut, tolerance)

    # display results
    print(f"a_1 : {a_mut_1:.8f}, b_1 : {b_mut_1:.8f}, c : {c_mut:.8f}")
    print(f"a_2 : {a_mut_2:.8f}, b_2 : {b_mut_2:.8f}, c : {c_mut:.8f}")

# Run Main function
main()

# Ankit Roy
# 22nd February, 2023
