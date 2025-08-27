#!/usr/bin/env python
# coding: utf-8

# The goal of this script is to plot the evolution of a two transfer functions that are interesting
# to understand CMB physics, the metric fluctuation phi and the effective temperature Theta_0+psi,
# as a function of both k and tau.


# import necessary modules
# uncomment to get plots displayed in notebook
#get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
import numpy as np
from classy import Class
from scipy.interpolate import interp1d
import math


#############################################
#
# User settings controlling the figure aspect
#
z_max_pk = 46000       # highest redshift involved
k_per_decade = 400     # number of k values, controls final resolution
k_min_tau0 = 40.       # this value controls the minimum k value in the figure (it is k_min * tau0)
P_k_max_inv_Mpc =1.0   # this value is directly the maximum k value in the figure in Mpc
tau_num_early = 2000   # number of conformal time values before recombination, controls final resolution
tau_num_late = 200     # number of conformal time values after recombination, controls final resolution
tau_ini = 10.          # first value of conformal time in Mpc
tau_label_Hubble = 20. # value of time at which we want to place the label on Hubble crossing
tau_label_ks = 40.     # value of time at which we want to place the label on sound horizon crossing
tau_label_kd = 230.    # value of time at which we want to place the label on damping scale crossing
#
# Cosmological parameters and other CLASS parameters
#
common_settings = {# LambdaCDM parameters (Planck 18 + lensing + BAO bestfit)
                   'omega_b':2.255065e-02,
                   'omega_cdm':1.193524e-01,
                   'H0':6.776953e+01,
                   'A_s':2.123257e-09,
                   'n_s':9.686025e-01,
                   'z_reio':8.227371e+00,
                   # As an output, we need the transfer functions only
                   'output':'mTk',
                   # Take fixed value for primordial Helium (instead of automatic BBN adjustment)
                   #'YHe':0.246,
                   # other output and precision parameters
                   'z_max_pk':z_max_pk,
                   'k_per_decade_for_pk':k_per_decade,
                   'k_per_decade_for_bao':k_per_decade,
                   'k_min_tau0':k_min_tau0, # this value controls the minimum k value in the figure
                   'perturbations_sampling_stepsize':'0.05',
                   'P_k_max_1/Mpc':P_k_max_inv_Mpc,
                   # At some point we will compute the evolution of the diffusion damping scale.
                   'compute damping scale':'yes',
                   # We want to display all transfer functions in the Newtonian gauge.
                   'gauge':'newtonian',
                   # We have set the gauge used internally in the code to be the Newtonian one.
                   # This line says that we want the source/transfer functions to be given by the code in this gauge,
                   # rather than being trasformed to gauge-invariant variables.
                   'matter_source_in_current_gauge':'yes'}

# create an instance of CLASS
M = Class()
# set input parameters
M.set(common_settings)


# Define the array used to sample conformal time in the final plot.
#
# For this, we first need to know the conformal time at recombination and today
# These values are extracted as derived parameters.
times = M.get_current_derived_parameters(['tau_rec','conformal_age'])
tau_rec=times['tau_rec']
tau_0 = times['conformal_age']
#
# define time sampling before recombination
tau1 = np.logspace(math.log10(tau_ini),math.log10(tau_rec),tau_num_early)
# define time sampling after recombination
tau2 = np.logspace(math.log10(tau_rec),math.log10(tau_0),tau_num_late)[1:]
tau2[-1] *= 0.999 # this tiny shift avoids interpolation errors
tau = np.concatenate((tau1,tau2))
tau_num = len(tau)


# use table of background and thermodynamics quantitites to define some functions
# returning some characteristic scales
# (of Hubble crossing, sound horizon crossing, etc.) at different times
#
background = M.get_background() # load background table
#print background.viewkeys()
thermodynamics = M.get_thermodynamics() # load thermodynamics table
#print thermodynamics.viewkeys()
#
background_tau = background['conf. time [Mpc]'] # read conformal times in background table
background_z = background['z'] # read redshift
background_aH = 2.*math.pi*background['H [1/Mpc]']/(1.+background['z'])/M.h() # read 2pi * aH in [h/Mpc]
background_ks = 2.*math.pi/background['comov.snd.hrz.']/M.h() # read 2pi/(comoving sound horizon) in [h/Mpc]
background_rho_m_over_r =\
    (background['(.)rho_b']+background['(.)rho_cdm'])\
    /(background['(.)rho_g']+background['(.)rho_ur']) # read rho_r / rho_m (to find time of equality)
background_rho_l_over_m =\
    background['(.)rho_lambda']\
    /(background['(.)rho_b']+background['(.)rho_cdm']) # read rho_m / rho_lambda (to find time of equality)
thermodynamics_tau = thermodynamics['conf. time [Mpc]'] # read confromal times in thermodynamics table
thermodynamics_kd = 2.*math.pi/thermodynamics['r_d']/M.h() # read 2pi(comoving diffusion scale) in [h/Mpc]
#
# define a bunch of interpolation functions based on previous quantities
#
background_z_at_tau = interp1d(background_tau,background_z)
background_aH_at_tau = interp1d(background_tau,background_aH)
background_ks_at_tau = interp1d(background_tau,background_ks)
background_tau_at_mr = interp1d(background_rho_m_over_r,background_tau)
background_tau_at_lm = interp1d(background_rho_l_over_m,background_tau)
thermodynamics_kd_at_tau = interp1d(thermodynamics_tau, thermodynamics_kd)
#
# infer arrays of characteristic quantitites calculated at values of conformal time in tau array
#
aH = background_aH_at_tau(tau)
ks = background_ks_at_tau(tau)
kd = thermodynamics_kd_at_tau(tau)
#
# infer times of R/M and M/Lambda equalities
#
tau_eq = background_tau_at_mr(1.)
tau_lambda = background_tau_at_lm(1.)
#
# check and inform user whether intial arbitrary choice of z_max_pk was OK
max_z_needed = background_z_at_tau(tau[0])
if max_z_needed > z_max_pk:
    print ('you must increase the value of z_max_pk to at least ',max_z_needed)
    () + 1  # this strange line is just a trick to stop the script execution there
else:
    print ('in a next run with the same values of tau, you may decrease z_max_pk from ',z_max_pk,' to ',max_z_needed)


# get transfer functions at each time and build arrays Theta0(tau,k) and phi(tau,k)
#
for i in range(tau_num):
    one_time = M.get_transfer(background_z_at_tau(tau[i])) # transfer functions at each time tau
    if i ==0:   # if this is the first time in the loop: create the arrays (k, Theta0, phi)
        k = one_time['k (h/Mpc)']
        k_num = len(k)
        Theta0 = np.zeros((tau_num,k_num))
        phi = np.zeros((tau_num,k_num))
    Theta0[i,:] = 0.25*one_time['d_g'][:]
    phi[i,:] = one_time['phi'][:]
#
# find the global extra of Theta0(tau,k) and phi(tau,k), used to define color code later
#
Theta_amp = max(Theta0.max(),-Theta0.min())
phi_amp = max(phi.max(),-phi.min())
#
# reshaping of (k,tau) necessary to call the function 'pcolormesh'
#
K,T = np.meshgrid(k,tau)
#
# inform user of the size of the grids (related to the figure resolution)
#
print ('grid size:',len(k),len(tau),Theta0.shape)


# start plotting
fig = plt.figure(figsize=(18,8))
#
# plot Theta0(k,tau)
#
ax_Theta = fig.add_subplot(121)
print ('> Plotting Theta_0')
fig_Theta = ax_Theta.pcolormesh(K,T,Theta0,cmap='coolwarm',vmin=-Theta_amp,vmax=Theta_amp,shading='auto')
print ('> Done')
#
# plot lines (characteristic times and scales)
#
ax_Theta.axhline(y=tau_rec,color='k',linestyle='-')
ax_Theta.axhline(y=tau_eq,color='k',linestyle='-')
ax_Theta.axhline(y=tau_lambda,color='k',linestyle='-')
ax_Theta.plot(aH,tau,'r-',linewidth=2)
ax_Theta.plot(ks,tau,color='#FFFF33',linestyle='-',linewidth=2)
ax_Theta.plot(kd,tau,'b-',linewidth=2)
#
# dealing with labels
#
ax_Theta.set_title(r'$\Theta_0$')
ax_Theta.text(1.5*k[0],0.9*tau_rec,r'$\mathrm{rec.}$')
ax_Theta.text(1.5*k[0],0.9*tau_eq,r'$\mathrm{R/M} \,\, \mathrm{eq.}$')
ax_Theta.text(1.5*k[0],0.9*tau_lambda,r'$\mathrm{M/L} \,\, \mathrm{eq.}$')
ax_Theta.annotate(r'$\mathrm{Hubble} \,\, \mathrm{cross.}$',
                  xy=(background_aH_at_tau(tau_label_Hubble),tau_label_Hubble),
                  xytext=(0.1*background_aH_at_tau(tau_label_Hubble),0.8*tau_label_Hubble),
                  arrowprops=dict(facecolor='black', shrink=0.05, width=1, headlength=5, headwidth=5))
ax_Theta.annotate(r'$\mathrm{sound} \,\, \mathrm{horizon} \,\, \mathrm{cross.}$',
                  xy=(background_ks_at_tau(tau_label_ks),tau_label_ks),
                  xytext=(0.07*background_aH_at_tau(tau_label_ks),0.8*tau_label_ks),
                  arrowprops=dict(facecolor='black', shrink=0.05, width=1, headlength=5, headwidth=5))
ax_Theta.annotate(r'$\mathrm{damping} \,\, \mathrm{scale} \,\, \mathrm{cross.}$',
                  xy=(thermodynamics_kd_at_tau(tau_label_kd),tau_label_kd),
                  xytext=(0.2*thermodynamics_kd_at_tau(tau_label_kd),2.0*tau_label_kd),
                  arrowprops=dict(facecolor='black', shrink=0.05, width=1, headlength=5, headwidth=5))
#
# dealing with axes
#
ax_Theta.set_xlim(k[0],k[-1])
ax_Theta.set_xscale('log')
ax_Theta.set_yscale('log')
ax_Theta.set_xlabel(r'$k  \,\,\, \mathrm{[h/Mpc]}$')
ax_Theta.set_ylabel(r'$\tau   \,\,\, \mathrm{[Mpc]}$')
ax_Theta.invert_yaxis()
#
# color legend
#
fig.colorbar(fig_Theta)
#
# plot phi(k,tau)
#
ax_phi = fig.add_subplot(122)
ax_phi.set_xlim(k[0],k[-1])
#ax_phi.pcolor(K,T,phi,cmap='coolwarm')
print ('> Plotting phi')
fig_phi = ax_phi.pcolormesh(K,T,phi,cmap='coolwarm',vmin=-0.,vmax=phi_amp,shading='auto')
print ('> Done')
#
# plot lines (characteristic times and scales)
#
ax_phi.axhline(y=tau_rec,color='k',linestyle='-')
ax_phi.axhline(y=tau_eq,color='k',linestyle='-')
ax_phi.axhline(y=tau_lambda,color='k',linestyle='-')
ax_phi.plot(aH,tau,'r-',linewidth=2)
ax_phi.plot(ks,tau,color='#FFFF33',linestyle='-',linewidth=2)
#
# dealing with labels
#
ax_phi.set_title(r'$\phi$')
ax_phi.text(1.5*k[0],0.9*tau_rec,r'$\mathrm{rec.}$')
ax_phi.text(1.5*k[0],0.9*tau_eq,r'$\mathrm{R/M} \,\, \mathrm{eq.}$')
ax_phi.text(1.5*k[0],0.9*tau_lambda,r'$\mathrm{M/L} \,\, \mathrm{eq.}$')
ax_phi.annotate(r'$\mathrm{Hubble} \,\, \mathrm{cross.}$',
                  xy=(background_aH_at_tau(tau_label_Hubble),tau_label_Hubble),
                  xytext=(0.1*background_aH_at_tau(tau_label_Hubble),0.8*tau_label_Hubble),
                  arrowprops=dict(facecolor='black', shrink=0.05, width=1, headlength=5, headwidth=5))
ax_phi.annotate(r'$\mathrm{sound} \,\, \mathrm{horizon} \,\, \mathrm{cross.}$',
                  xy=(background_ks_at_tau(tau_label_ks),tau_label_ks),
                  xytext=(0.07*background_aH_at_tau(tau_label_ks),0.8*tau_label_ks),
                  arrowprops=dict(facecolor='black', shrink=0.05, width=1, headlength=5, headwidth=5))
#
# dealing with axes
#
ax_phi.set_xscale('log')
ax_phi.set_yscale('log')
ax_phi.set_xlabel(r'$k \,\,\, \mathrm{[h/Mpc]}$')
ax_phi.set_ylabel(r'$\tau \,\,\, \mathrm{[Mpc]}$')
ax_phi.invert_yaxis()
#
# color legend
#
fig.colorbar(fig_phi)
#
# produce the PDF
#
#plt.show()
plt.savefig('many_times.png',dpi=300)
