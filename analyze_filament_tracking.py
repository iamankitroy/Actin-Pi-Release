#!/Users/roy/anaconda3/bin/python

"""

This script is used to analyze filament length data tracked with JFilament.
Please make sure that the JFilament output file is converted into a plottable CSV format.
The plottable CSV format should contain two columns: frame and length.
The script calculates the following parameters from the frame and length data:


Input: <filament-length-tracked-filename>.csv
Output: <filament-length-tracked-filename>_analysed.csv

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
import matplotlib.pyplot as plt

# Get input data
def dataIN(filename):
    data = pd.read_csv(filename)
    return data

# Normalize time and monomers
def normalizeProfile(data, peakTime, peakMonomers):

    # shift time axis
    data['norm_time'] = data.time - peakTime
    # shift monomer count axis
    data['norm_monomers'] = data.monomers - peakMonomers

    return data

# Calculate slope from a linear fit
def local_LinFit(x, y):
    m, b = np.polyfit(x, y, 1)
    return m

# Calculate instant polymerization velocity
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

# Write output data file
def dataOUT(data, filename):
    # output file name
    outname = '.'.join(filename.split('.')[:-1])
    outname = f"{outname}_analyzed.csv"
    
    data.to_csv(outname, index=False)

# Main function
def main():
    filename = sys.argv[1]
    label = filename.split('_')[0]

    time_interval = 5               # s
    pixel_size = 178                # in nm
    helical_rise = 2.75             # in nm

    elongation_data = dataIN(filename)

    elongation_data['time'] = elongation_data['frame'] * time_interval              # reaction time in s
    elongation_data['length_nm'] = elongation_data['length'] * pixel_size           # filament length in nm
    elongation_data['monomers'] = elongation_data['length_nm']/helical_rise         # number of monomers in filament

    # time with maximum filament length
    peakTime = int(elongation_data.loc[elongation_data.monomers == elongation_data.monomers.max(), 'time'])
    # maximum number of monomers in filament
    peakMonomers = float(elongation_data.monomers.max())

    # normalize elongation data to align filament peaks
    elongation_data = normalizeProfile(elongation_data, peakTime, peakMonomers)

    # label data
    elongation_data['label'] = label
    elongation_data['filename'] = filename
    
    # calculate instant growth velocities with entire
    velocities = calc_instantVelocity(elongation_data.time, elongation_data.monomers, 25)
    elongation_data['velocity'] = elongation_data['time'].map(velocities)

    dataOUT(elongation_data, filename)

# Run main
main()

# Ankit Roy
# 8th December, 2022
