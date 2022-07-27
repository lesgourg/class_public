#!/usr/bin/python

# This is an example for a nice way to calculate all primordial spectra
# and the background GW energy density in one python file.
import sys
import numpy as np
from scipy import constants


def read_parameter(argv=[]):
    par = {}
    # Standart values
    par['k_min']   = 1e-6
    par['k_max']   = 10
    par['k_size']  = 1000
    par['f_min']   = 0.5e-3
    par['f_max']   = 2e2
    par['f_size']  = 1000
    
    par['k_pivot'] = 0.05 # 1/Mpc
    par['A_s']     = 2.215e-09
    par['n_s']     = 0.9624
    par['r']       = 0.07
    par['n_t']     = -0.1
    par['A_gwi']   = 1e-10
    par['n_gwi']   = 0.
    par['c_s_gwi'] = 0.

    par['f_pivot']   = 1. # Hz
    par['Omega_gwb'] = 1e-5
    par['n_gwb']     = 0.
    
    # Read input parameter
    try:
        par['A_s']     = float(argv[0])
        par['n_s']     = float(argv[1])
        par['r']       = float(argv[2])
        par['n_t']     = float(argv[3])
        par['A_gwi']   = float(argv[4])
        par['n_gwi']   = float(argv[5])
        par['c_s_gwi'] = float(argv[6])

        par['Omega_gwb'] = float(argv[7])
        par['n_gwb']     = float(argv[8])
    except IndexError:
        pass
    except ValueError:
        raise ValueError("It seems some of the arguments are not correctly formatted. "+
                         "Remember that they must be floating point numbers.")
    
    return par

def init_ks(par):
    """Return the ks array in 1/Mpc"""
    return np.geomspace(par['k_min'], par['k_max'], par['k_size'])

def init_fs(par):
    """Return the fs array in Hz"""
    return np.geomspace(par['f_min'], par['f_max'], par['f_size'])


def f_to_k(f):
    """Convert frequency f to comoving wavenumber k"""
    return f * (2 * np.pi) / constants.c * (1e6 * constants.parsec)  # 1/Mpc

def k_to_f(k):
    """Convert comoving wavenumber k to frequency f"""
    return k / (2 * np.pi) * constants.c / (1e6 * constants.parsec)  # Hz


# Here we define the functions for the primordial spectra
def P_s(k, par):
    return par['A_s'] * (k/par['k_pivot'])**(par['n_s']-1.)

def P_t(k, par):
    return par['r']*par['A_s'] * (k/par['k_pivot'])**(par['n_t'])

def P_gwi(k, par):
    return par['A_gwi'] * (k/par['k_pivot'])**(par['n_gwi'])

def cross_s_gwi(k, par):
    return par['c_s_gwi']

# Background energy density of GWs
def Omega_GW(f, par):
    return par['Omega_gwb'] * (f/par['f_pivot'])**(par['n_gwb'])


# Print functions
def print_Pks(ks, par):
    print("# Dimensionless primordial spectrum, equal to [k^3/2pi^2] P(k)")
    print("# k [1/Mpc]                P_scalar(k)                P_gwi(k)                   ad x gwi                   P_tensor(k)")
    for k in ks:
        print("%.18e   %.18e   %.18e   %.18e   %.18e" % (k, P_s(k, par), P_gwi(k, par), cross_s_gwi(k, par), P_t(k, par)))

def print_Omega_GW(fs, par):
    print("# Dimensionless graviational wave background energy density Omega_GW(f)")
    print("# f [Hz]                   Omega_GW(f)")
    for f in fs:
        print("%.18e   %.18e" % (f, Omega_GW(f, par)))


def main(argv):
    if len(argv) == 1:
        raise Exception("Specify if you want to calculate the 'Pk' or 'OmGW'!")
        return 1

    if argv[1] == 'Pk':
        # Calaculate the P(k) spectrum
        # Initialize variables
        argv = argv[2:]
        par = read_parameter(argv)
        ks = init_ks(par)

        # Output spectrum
        print_Pks(ks, par)
    
    elif argv[1] == 'OmGW':
        # Calaculate the Omega_GW spectrum
        # Initialize variables
        argv = argv[2:]
        par = read_parameter(argv)
        fs = init_fs(par)

        # Output spectrum
        print_Omega_GW(fs, par)

    else:
        raise Exception("Can't interpret %s as input. You have to choose 'Pk' or 'OmGW'!" % argv[1])
        return 1

    return 0


if __name__ == '__main__':
    main(sys.argv)