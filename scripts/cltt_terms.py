#!/usr/bin/env python
# coding: utf-8

# The goal of this script is to plot the total CMB temperature spectrum C_l^TT,
# as well as its decomposition in different contributions:
# - T + SW: intrinsic temperature plus Sachs-Wolfe correction
# - early-ISW: early integrated Sachs-Wolfe
# - late-ISW: late integrated Sachs-Wolfe
# - Doppler contribution
# - total unlensed spectrum
# - total lensed spectrum


# import necessary modules
from classy import Class
from math import pi
#get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib
import matplotlib.pyplot as plt


####################################################
#
# Cosmological parameters and other CLASS parameters
#
####################################################
common_settings = {# LambdaCDM parameters (Planck 18 + lensing + BAO bestfit)
                   'omega_b':2.255065e-02,
                   'omega_cdm':1.193524e-01,
                   'H0':6.776953e+01,
                   'A_s':2.123257e-09,
                   'n_s':9.686025e-01,
                   'z_reio':8.227371e+00,
                   # output and precision parameters
                   'output':'tCl,pCl,lCl',
                   'lensing':'yes',
                   'l_max_scalars':5000}
###############
#
# Initiate a CLASS instance
#
M = Class()
#
###############
#
# call CLASS for the total Cl's and then for each contribution
#
###############
#
# total contribution
M.set(common_settings)
cl_tot = M.raw_cl(3000)   # raw_cl gives the unlensed spectrum
cl_lensed = M.lensed_cl(3000) # lensed_cl gives the lensed spectrum
#
# T+ISW contribution (unlensed)
M.empty() # reset input
M.set(common_settings) # new input
M.set({'temperature contributions':'tsw'})
cl_tsw = M.raw_cl(3000)
#
# early-ISW contribution (unlensed)
M.empty()
M.set(common_settings)
M.set({'temperature contributions':'eisw'})
cl_eisw = M.raw_cl(3000)
#
# late-ISW contribution (unlensed)
M.empty()
M.set(common_settings)
M.set({'temperature contributions':'lisw'})
cl_lisw = M.raw_cl(3000)
#
# Doppler contribution (unlensed)
M.empty()
M.set(common_settings)
M.set({'temperature contributions':'dop'})
cl_dop = M.raw_cl(3000)


#################
#
# start plotting
#
#################
#
plt.xlim([2,3000])
plt.xlabel(r"$\ell$")
plt.ylabel(r"$\ell (\ell+1) C_l^{TT} / 2 \pi \,\,\, [\times 10^{10}]$")
plt.grid()
#
ell = cl_tot['ell']
factor = 1.e10*ell*(ell+1.)/2./pi
plt.semilogx(ell,factor*cl_tsw['tt'],'c-',label=r'$\mathrm{T+SW}$')
plt.semilogx(ell,factor*cl_eisw['tt'],'r-',label=r'$\mathrm{early-ISW}$')
plt.semilogx(ell,factor*cl_lisw['tt'],'y-',label=r'$\mathrm{late-ISW}$')
plt.semilogx(ell,factor*cl_dop['tt'],'g-',label=r'$\mathrm{Doppler}$')
plt.semilogx(ell,factor*cl_tot['tt'],'r-',label=r'$\mathrm{total}$')
plt.semilogx(ell,factor*cl_lensed['tt'],'k-',label=r'$\mathrm{lensed}$')
#
plt.legend(loc='right',bbox_to_anchor=(1.4, 0.5))
plt.savefig('cltt_terms.pdf',bbox_inches='tight')
