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
    a = 0.18350639
    b = 0.01073622
    c = 0.15155544

    v_adp_depol = get_v_adp_depol(c)
    kr = get_kr(b)
    v_adp_pi_depol = get_adp_pi_depol(a, c)
    kr_be = get_kr_be(v_adp_depol, v_adp_pi_depol)

    print(f"{'v_ADP_depol':<15}: {v_adp_depol:.3f}")
    print(f"{'kr':<15}: {kr:.3f}")
    print(f"{'v_ADP-Pi_depol':<15}: {v_adp_pi_depol:.3f}")
    print(f"{'kr_BE':<15}: {kr_be:.3f}")
    print(f"{'kr_BE/kr':<15}: {kr_be/kr:.1f}")

# Run main
main()

# Ankit Roy
# 21st February, 2023
