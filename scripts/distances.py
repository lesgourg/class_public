#!/usr/bin/env python
# coding: utf-8

# The goal of this script is to compare the three cosmological distances:
# luminosity distance, angular diameter distance, comoving distance,
# to show how to extract background quantities from a CLASS run
# It creates a figure nearly identical to Figure 2.3 of Dodelsons's book 'Modern cosmology', first edition.


# import necessary modules
# uncomment to get plots displayed in notebook
#get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
import numpy as np
from classy import Class


# Create an instance of CLASS for the LambdaCDM model
LCDM = Class()
# Keep the default LCDM parameters but overwrite Omega_cdm, Omega_b; their sum is 0.3 such that Omega_Lambda=0.7
LCDM.set({'Omega_cdm':0.25,'Omega_b':0.05})


# Create an instance of CLASS for an Einstein-de Sitter model with the same parameters excepted
# Omega_cdm, Omega_b, now set such that their sum is 1, and thus Omega_Lambda=0
CDM = Class()
CDM.set({'Omega_cdm':0.95,'Omega_b':0.05})

# Just to cross-check that Omega_Lambda is negligible
# (but not exactly zero because we neglected radiation when we did our settings for Omega_cdm, Omega_b)
derived = CDM.get_current_derived_parameters(['Omega0_lambda'])
print (derived)
print ("Omega_Lambda =",derived['Omega0_lambda'])


# Get background quantities and check the list of keys in the output
baLCDM = LCDM.get_background()
baCDM = CDM.get_background()
baCDM.keys()


# Get H_0 in order to plot distances in this unit
fLCDM = LCDM.Hubble(0)
fCDM = CDM.Hubble(0)


# array of keys corresponding to the luminosity distance, comoving distance, and angular diameter distance
namelist = ['lum. dist.','comov. dist.','ang.diam.dist.']
colours = ['b','g','r']

# plot the three distances in the LambdaCDM model
for name in namelist:
    idx = namelist.index(name)
    plt.loglog(baLCDM['z'],fLCDM*baLCDM[name],colours[idx]+'-')


# plot the three distances in the Einstein-De Sitter model
for name in namelist:
    idx = namelist.index(name)
    plt.loglog(baCDM['z'],fCDM*baCDM[name],colours[idx]+'--')

plt.legend(namelist,loc='upper left')
plt.xlim([0.07, 10])
plt.ylim([0.08, 20])

plt.xlabel(r"$z$")
plt.ylabel(r"$\mathrm{Distance}\times H_0$")
plt.tight_layout()
plt.savefig('distances.pdf')
