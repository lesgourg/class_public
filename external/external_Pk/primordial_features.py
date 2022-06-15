#!/usr/bin/python

# This is an example for a nice way to calculate all primordial spectra
# and the background GW energy density in one python file.
import sys
import numpy as np


def read_parameter(argv=[]):
    par = {}
    # Standart values
    par['k_min']   = 1e-6
    par['k_max']   = 10
    par['k_size']  = 1000
    
    par['k_pivot'] = 0.05 # 1/Mpc
    par['A_s']     = 2.215e-09
    par['n_s']     = 0.9624
    par['alpha_s'] = 0.
    par['r']       = 0.07
    par['n_t']     = -0.1
    par['alpha_t'] = 0.

    par['A_lin']   = 0.
    par['w_lin']   = 0.
    par['phi_lin'] = 0.
    
    par['A_log']   = 0.
    par['w_log']   = 0.
    par['k_ref']   = par['k_pivot']
    par['phi_log'] = 0.
    
    # Read input parameter
    try:
        par['k_pivot'] = float(argv[0])
        par['A_s']     = float(argv[1])
        par['n_s']     = float(argv[2])
        # par['alpha_s'] = float(argv[3])
        par['r']       = float(argv[3])
        # par['n_t']     = float(argv[4])
        # par['alpha_t'] = float(argv[5])

        par['A_lin']   = float(argv[4])
        par['w_lin']   = float(argv[5])
        par['phi_lin'] = float(argv[6])

        par['A_log']   = float(argv[7])
        par['w_log']   = float(argv[8])
        par['k_ref']   = float(argv[9])
        par['phi_log'] = float(argv[6]) # = phi_lin
    except IndexError:
        pass
    except ValueError:
        raise ValueError("It seems some of the arguments are not correctly formatted. "+
                         "Remember that they must be floating point numbers.")
    
    return par

def init_ks(par):
    """Return the ks array in 1/Mpc"""
    return np.geomspace(par['k_min'], par['k_max'], par['k_size'])


# Here we define the functions for the primordial spectra
def P_s(k, par):
    x = np.log(k / par['k_pivot'])
    return par['A_s'] * np.exp((par['n_s']-1.) * x + 0.5 * par['alpha_s'] * x**2) \
            * (1. + par['A_lin'] * np.cos(par['w_lin']*k + par['phi_lin'])) \
            * (1. + par['A_log'] * np.cos(par['w_log']*np.log(k/par['k_ref']) + par['phi_log']))

def P_t(k, par):
    x = np.log(k / par['k_pivot'])
    return par['r']*par['A_s'] * np.exp(par['n_t'] * x + 0.5 * par['alpha_t'] * x**2) \
            * (1. + par['A_lin'] * np.cos(par['w_lin']*k + par['phi_lin'])) \
            * (1. + par['A_log'] * np.cos(par['w_log']*np.log(k/par['k_ref']) + par['phi_log']))


# Print functions
def print_Pks(ks, par):
    print("# Dimensionless primordial spectrum, equal to [k^3/2pi^2] P(k)")
    print("# k [1/Mpc]                P_scalar(k)                P_tensor(k)")
    for k in ks:
        print("%.18e   %.18e   %.18e" % (k, P_s(k, par), P_t(k, par)))


def main(argv):
    # Calaculate the P(k) spectrum
    # Initialize variables
    argv = argv[1:]
    par = read_parameter(argv)
    ks = init_ks(par)

    # Output spectrum
    print_Pks(ks, par)
    
    return 0


if __name__ == '__main__':
    main(sys.argv)