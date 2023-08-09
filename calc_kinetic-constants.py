#!/Users/roy/anaconda3/bin/python

"""

This script is used to calculate actin depolymerization kinetic constants.
Given parameters for an exponential decay model of 1/v_depol vs filament age, calculates the following:
v_ADP_depol		-> depolymerization velocity of ADP-actin
v_ADP-Pi_depol	-> depolymerization velocity of ADP-Pi-actin
kr				-> Pi release rate
kr_BE			-> Pi release rate from the barbed end
kr_BE/kr		-> kr_BE fold increase in comparison to kr

Exponential decay model is expressed as:
y = a*exp(-bx) + c
where:
x = filament age
y = 1/v_depol
a + c = 1/v_ADP-Pi_depol
b = kr
c = 1/v_ADP_depol

k_off_adp_pi is obtained from https://doi.org/10.1073/pnas.0702510104
kr_BE is obtained as described in https://doi.org/10.1371/journal.pbio.1001161

"""

__author__ = "Ankit Roy"
__copyright__ = "Copyright 2023, Bieling Lab, Max Planck Institute of Molecular Physiology"
__credits__ = ["Ankit Roy"]
__license__ = "GPLv3"
__maintainer__ = "Ankit Roy"
__email__ = "ankitroy.post@gmail.com"
__status__ = "Production"
__version__ = "1.0"

import sys

# actin-ADP depolymerisation rate
def get_v_adp_depol(c):
    return 1/c

# Internal Pi release rate
def get_kr(b):
    return b

# actin-ADP-Pi depolymerisation velocity
# Convolution of actin-ADP-Pi depolymerisation rate and barbed end Pi release rate
def get_adp_pi_depol(a, c):
    return 1/(a+c)

# Barbed end Pi release rate
def get_kr_be(v_adp_depol, v_adp_pi_depol):
    k_off_adp_pi = 0.2              # from Pollard's lab https://doi.org/10.1073/pnas.0702510104
    return (v_adp_depol * (k_off_adp_pi - v_adp_pi_depol)) / (v_adp_pi_depol - v_adp_depol)

# Main function
def main():
    a = 0.1696
    b = 0.0077
    c = 0.1515

    v_adp_depol = get_v_adp_depol(c)
    kr = get_kr(b)
    v_adp_pi_depol = get_adp_pi_depol(a, c)
    kr_be = get_kr_be(v_adp_depol, v_adp_pi_depol)

    print(f"{'v_ADP_depol':<15}: {v_adp_depol:.4f}")
    print(f"{'kr':<15}: {kr:.4f}")
    print(f"{'v_ADP-Pi_depol':<15}: {v_adp_pi_depol:.4f}")
    print(f"{'kr_BE':<15}: {kr_be:.4f}")
    print(f"{'kr_BE/kr':<15}: {kr_be/kr:.1f}")

# Run main
main()

# Ankit Roy
# 21st February, 2023
# 9th August, 2023
#		-->	Updated docstring.
