#!/usr/bin/env python
# coding: utf-8

# The goal of this script is to show the impact of the dark matter annihilation parameter p_ann
# on the various observables: CMB spectrum of temperature (C_l^TT) and polarisation (C_l^EE),
# matter power spectrum P(k,z), and CMB lensing spectrum C_l^phiphi
# (the latter two are actually not be affected by p_ann)


# import necessary modules
# uncomment to get plots displayed in notebook
#get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
import numpy as np
from classy import Class
from scipy.optimize import fsolve


############################################
#
# Varying parameter (others fixed to default)
#
# With the input syntax of class <= 2.9 we used: annihilation = 1.e-5 m^3/s/Kg
# With the new syntax this is equivalent to DM_annihilation_efficiency = 1.11e-22 m^3/s/J
# (the ratio is a factor (c/[1 m/s])**2 = 9.e16)
#
var_name = 'DM_annihilation_efficiency'
var_array = np.linspace(0,1.11e-22,5)
var_num = len(var_array)
var_legend = r'$p_\mathrm{ann}$'
var_figname = 'pann'
#
#############################################
#
# Fixed settings
#
common_settings = {# LambdaCDM parameters (Planck 18 + lensing + BAO bestfit)
                   'omega_b':2.255065e-02,
                   'omega_cdm':1.193524e-01,
                   'H0':6.776953e+01,
                   'A_s':2.123257e-09,
                   'n_s':9.686025e-01,
                   'z_reio':8.227371e+00,
                   # output and precision parameters
                   'output':'tCl,pCl,lCl,mPk',
                   'lensing':'yes',
                   'P_k_max_1/Mpc':3.0,
                   'l_switch_limber':9
                   }
#
# array of k values in 1/Mpc
#
kvec = np.geomspace(1e-4,3,num=1000)
legarray = []
twopi = 2.*np.pi
#
# Create figures
#
fig_Pk, ax_Pk = plt.subplots()
fig_TT, ax_TT = plt.subplots()
fig_EE, ax_EE = plt.subplots()
fig_PP, ax_PP = plt.subplots()
#
# create an instance of CLASS
M = Class()
#
# loop over varying parameter values
#
for i,var in enumerate(var_array):
    #
    print (' * Compute with %s=%e'%(var_name,var))
    #
    # deal with colors and legends
    #
    if i == 0:
        var_color = 'k'
        var_alpha = 1.
        legarray.append(r'ref. $\Lambda CDM$')
    else:
        var_color = 'r'
        var_alpha = 1.*i/(var_num-1.)
    if i == var_num-1:
        legarray.append(var_legend)
    #
    # fix input parameters in CLASS
    #
    M.set(common_settings)
    M.set({var_name:var})
    #
    # get Cls
    #
    clM = M.lensed_cl(2500)
    ll = clM['ell'][2:]
    clTT = clM['tt'][2:]
    clEE = clM['ee'][2:]
    clPP = clM['pp'][2:]
    #
    # get P(k) for common k values (k in 1/Mpc, P(k) in Mpc^3)
    #
    pkM = M.get_pk_all(kvec,0.)
    #
    # get reduced Hubble h for unit conversion
    h = M.h()
    #
    # plot P(k) with k in h/Mpc and P(k) in (Mpc/h)^3
    #
    ax_Pk.loglog(kvec/h,np.array(pkM)*h*h*h,color=var_color,alpha=var_alpha,linestyle='-')
    #
    # plot C_l^TT
    #
    ax_TT.semilogx(ll,clTT*ll*(ll+1)/twopi,color=var_color,alpha=var_alpha,linestyle='-')
    #
    # plot C_l^EE
    #
    ax_EE.loglog(ll,clEE*ll*(ll+1)/twopi,color=var_color,alpha=var_alpha,linestyle='-')
    #
    # plot C_l^phiphi
    #
    ax_PP.loglog(ll,clPP*ll*(ll+1)*ll*(ll+1)/twopi,color=var_color,alpha=var_alpha,linestyle='-')
    #
    # reset CLASS
    #
    M.empty()
#
# output of P(k) figure
#
ax_Pk.set_xlim([1e-4,3.])
ax_Pk.set_xlabel(r'$k \,\,\,\, [h/\mathrm{Mpc}]$')
ax_Pk.set_ylabel(r'$P(k) \,\,\,\, [\mathrm{Mpc}/h]^3$')
ax_Pk.legend(legarray)
fig_Pk.tight_layout()
fig_Pk.savefig('varying_%s_Pk.pdf' % var_figname)
#
# output of C_l^TT figure
#
ax_TT.set_xlim([2,2500])
ax_TT.set_xlabel(r'$\ell$')
ax_TT.set_ylabel(r'$[\ell(\ell+1)/2\pi]  C_\ell^\mathrm{TT}$')
ax_TT.legend(legarray)
fig_TT.tight_layout()
fig_TT.savefig('varying_%s_cltt.pdf' % var_figname)
#
# output of C_l^EE figure
#
ax_EE.set_xlim([2,2500])
ax_EE.set_xlabel(r'$\ell$')
ax_EE.set_ylabel(r'$[\ell(\ell+1)/2\pi]  C_\ell^\mathrm{EE}$')
ax_EE.legend(legarray)
fig_EE.tight_layout()
fig_EE.savefig('varying_%s_clee.pdf' % var_figname)
#
# output of C_l^pp figure
#
ax_PP.set_xlim([10,2500])
ax_PP.set_xlabel(r'$\ell$')
ax_PP.set_ylabel(r'$[\ell^2(\ell+1)^2/2\pi]  C_\ell^\mathrm{\phi \phi}$')
ax_PP.legend(legarray)
fig_PP.tight_layout()
fig_PP.savefig('varying_%s_clpp.pdf' % var_figname)
