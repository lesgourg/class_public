#!/usr/bin/env python
# coding: utf-8

# The goal of this script is to plot the non-linear matter power spectrum computed with HMcode 2020.
# We also show for comparison the linear spectrum, as well as the no-wiggle power spectrum computed in two ways.


# import classy module
from classy import Class
# uncomment to get plots displayed in notebook
#get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
import numpy as np


# create an instance of Class
LambdaCDM = Class()
# pass input parameters
LambdaCDM.set_baseline('p18') # Start from Planck 2018 baseline
LambdaCDM.set({'output':'tCl,pCl,lCl,mPk','analytic_nowiggle':'yes','numerical_nowiggle':'yes','lensing':'yes'})
LambdaCDM.set({'P_k_max_h/Mpc':5.0,'z_max_pk':1.1})
LambdaCDM.set({'non_linear':'HMcode','hmcode_version':'2020'})


# create an array of k values in h/Mpc
kk = np.geomspace(1e-4,3,num=1000)
# execute CLASS such that we can extract the Hubble rate from the class instance
LambdaCDM.compute()
# get reduced Hubble for conversions to 1/Mpc
h = LambdaCDM.h()


# get P(k) at redhsift z
z=0.
# get linear P(k,z) in Mpc^3 and convert it to (Mpc/h)^3
Pk_lin = LambdaCDM.get_pk_all(kk*h,z,nonlinear=False)*h**3
# get non-linear P(k,z) in Mpc^3 and convert it to (Mpc/h)^3
Pk = LambdaCDM.get_pk_all(kk*h,z,nonlinear=True)*h**3
# get the linear no-wiggle spectrum from a numerical method (smoothing of the full spectrum)
Pk_nw = np.array([LambdaCDM.pk_numerical_nw(k*h,z)*h**3 for k in kk]) # note: this function will soon be upated to the get_pk_all framework
# get the linear no-wiggle spectrum from an analytic approximation (Eisenstein & Hu)
Pk_an_nw = np.array([LambdaCDM.pk_analytic_nw(k*h)*h**3 for k in kk]) # note: this function will soon be upated to the get_pk_all framework


# plot P(k)
#plt.figure(1)
plt.xscale('log');plt.yscale('log')
#plt.xlim(kk[0],kk[-1])
plt.xlim(1.e-3,0.5)
plt.ylim(200,3e4)
plt.xlabel(r'$k \,\,\,\, [h/\mathrm{Mpc}]$')
plt.ylabel(r'$P(k) \,\,\,\, [\mathrm{Mpc}/h]^3$')
plt.plot(kk,Pk_nw,'k-',label='linear, no-wiggle (smoothing)')
plt.plot(kk,Pk_an_nw,'r-',label='linear, no-wiggle (EH)')
plt.plot(kk,Pk_lin,'g-',label='linear')
plt.plot(kk,Pk,'b-',label='non-linear')
plt.legend()
plt.savefig('test_hmcode.pdf',bbox_inches='tight')
