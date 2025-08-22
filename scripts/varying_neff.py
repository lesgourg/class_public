#!/usr/bin/env python
# coding: utf-8

# The goal of this script is to show the impact of the effective neutrino number N_eff
# on the various observables: CMB spectrum of temperature (C_l^TT) and polarisation (C_l^EE),
# matter power spectrum P(k,z), and CMB lensing spectrum C_l^phiphi.
#
# We choose to vary at the same time (Neff, h, omega_cdm) in order to illustrate the partial
# degeneracy bwteeen these parameters: we vary them in such way to keep the ratio of
# radiation-to-matter and matter-to-Lambda energy fixed (that is, fixed z_eq and z_Lambda).


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
var_name = 'Neff'
var_array = np.linspace(3.044,5.044,5)
var_num = len(var_array)
var_legend = r'$N_\mathrm{eff}$'
var_figname = 'neff'
#
# Constraints to be matched in order to illustrate the partial degeneracy between N
#
# As explained in the "Neutrino cosmology" book, CUP, Lesgourgues et al., section 5.3, the goal is to vary
# - omega_cdm by a factor alpha = (1 + coeff*Neff)/(1 + coeff*3.046)
# - h by a factor sqrt*(alpha)
# in order to keep a fixed z_equality(R/M) and z_equality(M/Lambda)
#
# LambdaCDM parameters (Planck 18 + lensing + BAO bestfit)
omega_b = 2.255065e-02
omega_cdm_standard = 1.193524e-01
h_standard = 0.6776953
#
# coefficient such that omega_r = omega_gamma (1 + coeff*Neff),
# i.e. such that omega_ur = omega_gamma * coeff * Neff:
# coeff = omega_ur/omega_gamma/Neff_standard
# We could extract omega_ur and omega_gamma on-the-fly within th script,
# but for simplicity we did a preliminary interactive run with background_verbose=2
# and we copied the values given in the budget output.
#
coeff = 1.70961e-05/2.47298e-05/3.044
print ("coeff=",coeff)
#
#############################################
#
# Fixed settings
#
common_settings = {# fixed LambdaCDM parameters (Planck 18 + lensing + BAO bestfit)
                   'omega_b':omega_b,
                   'A_s':2.100549e-09,
                   'n_s':0.9660499,
                   'tau_reio':0.05430842,
                   # output and precision parameters
                   'output':'tCl,pCl,lCl,mPk',
                   'lensing':'yes',
                   'P_k_max_1/Mpc':3.0}
#
##############################################
#
# loop over varying parameter values
#
M = {}
#
for i, Neff in enumerate(var_array):
    #
    # rescale omega_cdm and h
    #
    alpha = (1.+coeff*Neff)/(1.+coeff*3.044)
    omega_cdm = (omega_b + omega_cdm_standard)*alpha - omega_b
    h = h_standard*np.sqrt(alpha)
    print (' * Compute with %s=%e, %s=%e, %s=%e'%('Neff',Neff,'omega_cdm',omega_cdm,'h',h))
    #
    # create and instance of CLASS
    M[i] = Class()
    # set input parameters
    M[i].set(common_settings)
    #M[i].set(better_precision_settings) # you could put also better precision settings if you wanted
    M[i].set({'Neff':Neff})
    M[i].set({'omega_cdm':omega_cdm})
    M[i].set({'h':h})


#############################################
#
# extract spectra and plot them
#
#############################################
#
# create array of kvec in h/Mpc
kvec = np.logspace(-4,np.log10(3),1000)
twopi = 2.*np.pi
#
# Create figures
#
fig_Pk, ax_Pk = plt.subplots()
fig_TT, ax_TT = plt.subplots()
#
# loop over varying parameter values
#
ll = {}
clM = {}
clTT = {}
pkM = {}
legarray = []

for i, Neff in enumerate(var_array):
    #
    # compute rescaling factor
    alpha = (1.+coeff*Neff)/(1.+coeff*3.044)
    # compute rescaled h for unit conversion
    h = 0.67810*np.sqrt(alpha)
    #
    # deal with colors and legends
    #
    if i == 0:
        var_color = 'k'
        var_alpha = 1.
    else:
        var_color = plt.cm.Reds(0.8*i/(var_num-1))
    #
    # get lensed Cls
    #
    clM[i] = M[i].lensed_cl(2500)
    # extract multipole values
    ll[i] = clM[i]['ell'][2:]
    # extract C_l^TT values
    clTT[i] = clM[i]['tt'][2:]
    #
    # get P(k), converted to units of (Mpc/h)^3
    # (note that kvec*h is the array of k values in 1/Mpc)
    #
    pkM[i] = M[i].get_pk_all(kvec*h,0.)*h**3
    #
    # plot P(k)
    #
    if i == 0:
        ax_Pk.semilogx(kvec,np.array(pkM[i])/np.array(pkM[0]),
                       color=var_color,#alpha=var_alpha,
                       linestyle='-')
    else:
        ax_Pk.semilogx(kvec,np.array(pkM[i])/np.array(pkM[0]),
                       color=var_color,#alpha=var_alpha,
                       linestyle='-',
                       label=r'$\Delta N_\mathrm{eff}=%g$'%(Neff-3.044))
    #
    # plot C_l^TT
    #
    if i == 0:
        ax_TT.semilogx(ll[i],clTT[i]/clTT[0],
                       color=var_color,alpha=var_alpha,linestyle='-')
    else:
        ax_TT.semilogx(ll[i],clTT[i]/clTT[0],
                       color=var_color,alpha=var_alpha,linestyle='-',
                       label=r'$\Delta N_\mathrm{eff}=%g$'%(Neff-3.044))
#
# output of P(k) figure
#
ax_Pk.set_xlim([1.e-3,3.])
ax_Pk.set_ylim([0.98,1.20])
ax_Pk.set_xlabel(r'$k \,\,\,\, [h^{-1}\mathrm{Mpc}]$')
ax_Pk.set_ylabel(r'$P(k)/P(k)[N_\mathrm{eff}=3.046]$')
ax_Pk.legend(loc='upper left')
fig_Pk.tight_layout()
fig_Pk.savefig('ratio-%s-Pk.pdf' % var_figname)
#
# output of C_l^TT figure
#
ax_TT.set_xlim([2,2500])
ax_TT.set_ylim([0.850,1.005])
ax_TT.set_xlabel(r'$\mathrm{Multipole} \,\,\,\,  \ell$')
ax_TT.set_ylabel(r'$C_\ell^\mathrm{TT}/C_\ell^\mathrm{TT}(N_\mathrm{eff}=3.046)$')
ax_TT.legend(loc='lower left')
fig_TT.tight_layout()
fig_TT.savefig('ratio-%s-cltt.pdf' % var_figname)
