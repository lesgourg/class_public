#!/usr/bin/env python
# coding: utf-8

# The goal of this script is to plot the evolution of the free electron fraction as a function of
# conformal time, to show how to extract thermodynamical quantities from a CLASS run


# import necessary modules
# uncomment to get plots displayed in notebook#%matplotlib inline
import matplotlib.pyplot as plt
import numpy as np
from classy import Class


# Create a CLASS instance
M = Class()

# We could customize our input parameters but here we choose to stick to the default LambdaCDM parameters
# No need to pass the 'output' field since we only need CLASS to compute the background and thermodynamics evolution 

# get derived parameters computed by CLASS for this cosmology: here we get
# the conformal time at recombination (tau_rec), today (conformal_age), and at reionization (conf_time_reio) 
derived = M.get_current_derived_parameters(['tau_rec','conformal_age','conf_time_reio'])

# extract all thermodynamics quantities at each redshift, sorted by growing redhsift, and print the list of available keys
thermo = M.get_thermodynamics()
print (thermo.keys())

# within the thermo dictionary, extract the array of conformal time tau
tau = thermo['conf. time [Mpc]']
# within the thermo dictionary, extract the array of visibility function g
g = thermo['g [Mpc^-1]']
# to make the reionisation peak visible, rescale g by 100 only for late times 
# (here we do it for the first 500 values of redshift, corresponding to the 500 largest values of tau)
g[:500] *= 100


# plot g(tau)    
plt.xlim([1.e2,derived['conformal_age']])
plt.xlabel(r'$\tau \,\,\, \mathrm{[Mpc]}$')
plt.ylabel(r'$\mathrm{visibility} \,\,\, g \,\,\, [\mathrm{Mpc}^{-1}]$')
plt.axvline(x=derived['tau_rec'],color='k')
plt.axvline(x=derived['conf_time_reio'],color='k')
#
plt.semilogx(tau,g,'r',label=r'$\psi$')
plt.savefig('thermo.pdf',bbox_inches='tight')

