#!/usr/bin/env python
# coding: utf-8

# The goal of this script is to plot the various CMB temperature spectra:
# - temperature spectrum C_l^TT from scalars and tensors
# - polarisation spectrum C_l^EE from scalars and tensors
# - polarisation spectrum C_l^BB from lesned scalars and tensors


# import necessary modules
from classy import Class
from math import pi
#get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib
import matplotlib.pyplot as plt


#####################################################
#
# Cosmological parameters and other CLASS parameters
#
#####################################################
common_settings = {# LambdaCDM parameters (Planck 18 + lensing + BAO bestfit)
                   'omega_b':2.255065e-02,
                   'omega_cdm':1.193524e-01,
                   'H0':6.776953e+01,
                   'A_s':2.123257e-09,
                   'z_reio':8.227371e+00}

# We compute scalars up to multipole l=3000
l_max_scalars = 3000
# Since tensors get negligible at large l, we compute them only up to l=600
l_max_tensors = 600

# Note that for l_max_tensors=600 we can keep default precision,
# while for accurate tensor predictions up to  l_max_tensors=3000
# we would need to import many high precision settings from the file cl_ref.pre


###############
#
# call CLASS : scalars only, with a given value of the scalar tilt n_s
#
###############
#
M = Class()
M.set(common_settings)
M.set({'output':'tCl,pCl','modes':'s','lensing':'no','n_s':0.9660499,
       'l_max_scalars':l_max_scalars})
cls = M.raw_cl(l_max_scalars)  # raw_cl gives the unlensed spectrum


###############
#
# call CLASS : tensors only, with a given value of the tensor tilt n_t
#              and of the tensor-to-scalar ratio r
#              (so A_t = r * A_s where A_s was passed before)
#
###############
#
M.empty() # reset input
M.set(common_settings) # new input
M.set({'output':'tCl,pCl','modes':'t','lensing':'no','r':0.1,'n_t':0,
       'l_max_tensors':l_max_tensors})
clt = M.raw_cl(l_max_tensors)  # raw_cl gives the unlensed spectrum


###############
#
# call CLASS : scalars + tensors (only in this case we can get the correct lensed C_l^BB)
#
###############
#
M.empty() # reset input parameters to default, before passing a new parameter set
M.set(common_settings)
M.set({'output':'tCl,pCl,lCl','modes':'s,t','lensing':'yes','n_s':0.9660499,'r':0.1,'n_t':0,
       'l_max_scalars':l_max_scalars,'l_max_tensors':l_max_tensors})
cl_tot = M.raw_cl(l_max_scalars)   # raw_cl gives the unlensed spectrum
cl_lensed = M.lensed_cl(l_max_scalars)  # lensed_cl gives the lensed spectrum


#################
#
# plotting
#
#################
#
plt.xlim([2,l_max_scalars])
plt.ylim([1.e-8,10])
plt.xlabel(r"$\ell$")
plt.ylabel(r"$\ell (\ell+1) C_l^{XY} / 2 \pi \,\,\, [\times 10^{10}]$")
plt.title(r"$r=0.1$")
plt.grid()
#
ell = cl_tot['ell']
ellt = clt['ell']
factor = 1.e10*ell*(ell+1.)/2./pi
factort = 1.e10*ellt*(ellt+1.)/2./pi
#
plt.loglog(ell,factor*cls['tt'],'r-',label=r'$\mathrm{TT(s)}$')
plt.loglog(ellt,factort*clt['tt'],'r:',label=r'$\mathrm{TT(t)}$')
plt.loglog(ell,factor*cls['ee'],'b-',label=r'$\mathrm{EE(s)}$')
plt.loglog(ellt,factort*clt['ee'],'b:',label=r'$\mathrm{EE(t)}$')
plt.loglog(ellt,factort*clt['bb'],'g:',label=r'$\mathrm{BB(t)}$')
plt.loglog(ell,factor*(cl_lensed['bb']-cl_tot['bb']),'g-',label=r'$\mathrm{BB(lensing)}$')
plt.legend(loc='right',bbox_to_anchor=(1.4, 0.5))
plt.savefig('cl_ST.pdf',bbox_inches='tight')
