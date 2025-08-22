#!/usr/bin/env python
# coding: utf-8

# The goal of this script is to compare various approximations for the
# growth factor D(a) and growth rate f(a)
# in a cosmology with Chevalier-Polarski-Linder (CPL) dynamical dark energy


#get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
import numpy as np
from classy import Class


# We want to use a few CPL fluid models, parameterized as usual through w0,wa
w0vec = [-0.7, -1.0, -1.3]
#wavec = [-0.2,0.0,0.2]
wavec = [0.0]
# If you want to test for just a cosmological constant, use the below
#w0vec = [-1.0]
#wavec = [0.0]


# Now, set cosmologies for each of the (w0,wa) combinations
cosmo = {}
for w0 in w0vec:
    for wa in wavec:
        # Give it a nice name:
        if w0==-1.0 and wa==0.0:
            M='LCDM'
        else:
            M = '('+str(w0)+','+str(wa)+')'
        # Then assign the correct parameters
        # If the model is LCDM, we could run it as a fluid with w0_fld = -1, wa_fld=0 ,
        # but it's probably nicer to do it as the usual cosmological constant
        # We also set the Gauge to be Newtonian, though this isn't critical
        cosmo[M] = Class()
        cosmo[M].set({'input_verbose':1,'background_verbose':1,'gauge' : 'Newtonian'})
        if M!='LCDM':
            cosmo[M].set({'Omega_Lambda':0.,'w0_fld':w0,'wa_fld':wa})


# To interpolate stuff
from scipy.interpolate import CubicSpline
# To compute analytical/numerical solutions
import scipy.special
import scipy.integrate


# Define hypergeometric approximation based on analytical arguments
def D_hypergeom(avec,csm):
    # From the background we want to get the Omega_DE(z=0), which is index -1
    bg = csm.get_background()
    # We also want the Omega_m today
    Om = csm.Omega0_m()
    # If there is a cosmological constant (lambda), take that as DE
    if '(.)rho_lambda' in bg:
        Ol = bg['(.)rho_lambda'][-1]/bg['(.)rho_crit'][-1]
    # Otherwise, assume there is a fluid, and take that as DE
    else:
        Ol = bg['(.)rho_fld'][-1]/bg['(.)rho_crit'][-1]

    # This is the analytically known result
    x = Ol/Om*avec**3
    D = avec*scipy.special.hyp2f1(1./3.,1,11./6.,-x)
    D_today = scipy.special.hyp2f1(1./3.,1,11./6.,-Ol/Om)
    return D/D_today

# Similarly to the growth factor D, the growth factor f also has an analytical result
def f_hypergeom(avec,csm):
    # Same as for the D_hypergeom
    bg = csm.get_background()
    Om = csm.Omega0_m()
    if '(.)rho_lambda' in bg:
        Ol = bg['(.)rho_lambda'][-1]/bg['(.)rho_crit'][-1]
    else:
        Ol = bg['(.)rho_fld'][-1]/bg['(.)rho_crit'][-1]

    # These are the analytically known formulae
    x = Ol/Om*avec**3
    D = avec*scipy.special.hyp2f1(1./3.,1,11./6.,-x)
    f = 1.-6./11.*x*avec/D*scipy.special.hyp2f1(4./3.,2,17./6.,-x)
    return f

# We can also try to compute results in a numerically approximate way
def D_integral2(avec,csm):
    # Again, we will need Omega_DE(z=0), computed as below
    # Additionally, we need w(z), which we use the CPL approximation w(z)=w0+wa*(1-a) for
    bg = csm.get_background()
    Om = csm.Omega0_m()
    # If there is a cosmological constant (lambda), take that as DE
    if '(.)rho_lambda' in bg:
        Ol = bg['(.)rho_lambda'][-1]/bg['(.)rho_crit'][-1]
        w0 = -1
        wa = 0.0
    # Otherwise, assume there is a fluid, and take that as DE
    else:
        Ol = bg['(.)rho_fld'][-1]/bg['(.)rho_crit'][-1]
        w0 = csm.pars['w0_fld']
        wa = csm.pars['wa_fld']
    # Now compute the growth function D by solving for each value of a the explicit integral equations
    D = np.zeros(avec.shape)
    # Define our approximation for the comoving Hubble parameter
    # We neglect radiation since we only integrat at low redshifts
    # We also don't include the H0 normalization since that would cancel out at the end anyways
    Hc = lambda a:a*np.sqrt(Om/a**3 + Ol*a**(-3*(1+w0+wa))*np.exp(-3.*(1.0-a)*wa) )
    # The integrand for computing how D(a) deviates from H(a) turns out to be I = integral 1/(aH)^3 da
    # The real equation is of course D(a) = 5/2 * Omega_m * H(a) * integral_0^a da/(a*H(a))**3
    Dintegrand2 = lambda a: Hc(a)**(-3)
    for idx, a in enumerate(avec):
        I = scipy.integrate.quad(Dintegrand2, 1e-4, a)
        D[idx] = Hc(a)/a*I[0]
    # Normalize by the value today -- this is where H0 would cancel out anyways
    D = D/scipy.integrate.quad(Dintegrand2,1e-4,1)[0]
    return D

# We can also try to compute results in a numerically approximate way
# This time, using a different integral approximation
def D_integral(avec,csm):
    # Again, we will need Omega_DE(z=0), computed as below
    bg = csm.get_background()
    Om = csm.Omega0_m()
    # Works only in LCDM where we have Lambda, a cosmological constant
    Ol = bg['(.)rho_lambda'][-1]/bg['(.)rho_crit'][-1]
    # This time, we take into account radiation
    Or = 1-Om-Ol
    # Our integrand now has radiation, compared to the above
    Hc = lambda a:np.sqrt(Om/a+Ol*a*a+Or/a/a)
    Dintegrand = lambda a:Hc(a)**(-3)
    D = np.zeros(avec.shape)
    for idx, a in enumerate(avec):
        # To avoid numerical issues, currently not needed:
        #if a<1e-4:
        #    continue
        # Now compute again the same integrals, this time using the different approxmation
        I = scipy.integrate.quad(Dintegrand,1e-15,a,args=())
        D[idx] = Hc(a)/a*I[0]
    # Again, we normalize by the value today
    D = D/scipy.integrate.quad(Dintegrand,1e-15,1,args=())[0]
    return D

# The approximation by Linder et al.
def D_linder(avec,csm):
    # This time, we once again need Omega_DE(z=0) and w(z)
    bg = csm.get_background()
    if '(.)rho_lambda' in bg:
        Ol = bg['(.)rho_lambda'][-1]/bg['(.)rho_crit'][-1]
        w0 = -1
        wa = 0.0
    else:
        Ol = bg['(.)rho_fld'][-1]/bg['(.)rho_crit'][-1]
        w0 = csm.pars['w0_fld']
        wa = csm.pars['wa_fld']

    # We also this time need Omega_m(a)
    # In principle it would be preferable to use
    #Om_of_a = (bg['(.)rho_cdm']+bg['(.)rho_b'])/bg['H [1/Mpc]']**2
    Om_of_a = csm.Om_m(bg['z'])
    a_bg = 1./(1.+bg['z'])

    # This is Linder's integrand approximation involving f(a) approximately Om(a)**gamma, and integrating that since dlnD(a)/dlna = f(a)
    gamma = 0.55+0.05*(w0+0.5*wa)
    integ = (Om_of_a**gamma-1.)/a_bg

    # For quickly computing the integrand at arbitrary values of a
    integ_interp = CubicSpline(a_bg,integ)

    # Now start to define a function to compute the desired integrals
    D = np.zeros(avec.shape)
    # In this case, we probably are fine with z=1/a-1=1/(1e-3)-1=999 as our maximum redshift
    # If you want to push it to the limit, you could also go all the way to the a_min of the background module
    # But honestly that's not advised
    #amin = min(a_bg)
    amin = 1e-3
    for idx, a in enumerate(avec):
        # Before amin we can approximate D(a) = a
        if a<amin:
            D[idx] = a
        # Afterwards we have to compute the integral
        else:
            I = scipy.integrate.quad(integ_interp,amin,a,args=())
            D[idx] = a*np.exp(I[0])
    # We are in a lucky case here: No re-normalization is required
    # This is an artifact of how the integrand is defined (with I=1 at a=1)
    #D = D/scipy.integrate.quad(Dintegrand,1e-15,1,args=())[0]
    return D

# The other integral based on Linder's approximation, this time ignoring all other components other than DE and matter
def D_linder2(avec,csm):
    # This time, we once again need Omega_DE(z=0) and w(z)
    # However, this time we also need rho_DE(z) and rho_M(z)
    bg = csm.get_background()
    if '(.)rho_lambda' in bg:
        Ol = bg['(.)rho_lambda'][-1]/bg['(.)rho_crit'][-1]
        w0 = -1
        wa = 0.0
        rho_de = bg['(.)rho_lambda']
    else:
        Ol = bg['(.)rho_fld'][-1]/bg['(.)rho_crit'][-1]
        w0 = csm.pars['w0_fld']
        wa = csm.pars['wa_fld']
        rho_de = bg['(.)rho_fld']

    rho_M = bg['(.)rho_cdm']+bg['(.)rho_b']

    # Instead of using the 'correct' formula below
    #Om_of_a = rho_M/bg['H [1/Mpc]']**2
    # We approximate assuming only DE + matter are present
    Om_of_a = rho_M/(rho_M+rho_de)

    # Otherwise, the idea is the same as in the Linder case above
    gamma = 0.55+0.05*(1+w0+0.5*wa)

    # We could do a fine integration again, but in this case
    # the results are very approximate anyway. So we can also
    # just do a trapezoidal integration over the supplied avec
    #a_bg = 1./(1.+bg['z'])
    a_bg = avec
    integ = (Om_of_a**gamma-1.)/a_bg

    # Do a trapezoidal integration in this case instead of an expensive full one
    # since we are anyways strongly approximating things
    D = np.zeros(avec.shape)
    for idx, a in enumerate(avec):
        if idx<2:
            I=0
        else:
            I = np.trapz(integ[:idx],x=avec[:idx])
        D[idx] = a*np.exp(I)
    # Again, we don't need to renormalize because of the way things are defined (with I=1 at a=1)
    #D = D/scipy.integrate.quad(Dintegrand,1e-15,1,args=())[0]
    return D/D[-1]



# This is just a convenience function for drawing a vertical line at a given redshift
# In particular, it converts a given z value into a variable of choice (like z, a, or tau) and then plots that as a vertical line
def draw_vertical_redshift(csm, theaxis, var='tau',z=99,ls='-.',label='$z=99$'):
    # Determine x value by converting z -> x = (z,a,or tau)
    if var=='z':
        xval = z
    elif var=='a':
        xval = 1./(z+1)
    elif var=='tau':
        bg = csm.get_background()
        f = CubicSpline(bg['z'],bg['conf. time [Mpc]'])
        xval = f(z)
    # Make a line at that x value with the desired label
    theaxis.axvline(xval,lw=1,ls=ls,color='k',label=label)


# Define figure properties
figwidth1 = 4.4 #=0.7*6.3
figwidth2 = 6.3
figwidth15 = 0.5*(figwidth1+figwidth2)
ratio = 8.3/11.7
figheight1 = figwidth1*ratio
figheight2 = figwidth2*ratio
figheight15 = figwidth15*ratio

# Define line properties
lw=2
fs=12
labelfs=16

# Create figure and axes
fig, (ax1, ax2) = plt.subplots(2,1,figsize=(1.2*figwidth1,figheight1/(3./5.)),sharex=True,
                              gridspec_kw = {'height_ratios':[3, 2]})

# Select the desired amin to be explored, currently just the lowest one is selected
if False:
    aminexp = -13
    amin = 10**aminexp
    ymin = 10**(aminexp/2.)
    ymax = 10**(-aminexp/2.)
elif False:
    aminexp = -7
    amin = 10**aminexp
    ymin = 10**(aminexp)
    ymax = 10**(-aminexp)
else:
    aminexp = -4
    amin = 10**aminexp
    ymin = 10**(aminexp-1)
    ymax = 10**(-aminexp+1)

# Get the background from LCDM
bg = cosmo['LCDM'].get_background()

# Get the default quantities in LCDM
a = 1./(bg['z']+1)
H = bg['H [1/Mpc]']
D = bg['gr.fac. D']
f = bg['gr.fac. f']

# Compare within LCDM the CLASS computation to the analytical result
ax1.loglog(a,D,lw=lw,label=r'$D_+^\mathrm{class}$')
ax1.loglog(a,D_hypergeom(a,cosmo['LCDM']),lw=lw,label=r'$D_+^\mathrm{analytic}$')

# Compare the results with known limits (plotted somewhere very different, in order not to disturb the other curves)
ax1.loglog(a,a*ymax,'k--',lw=lw,label=r'$\propto a$')
ax1.loglog(a,1./a*ymin,'k:',lw=lw,label=r'$\propto a^{-1}$')

# Also show the ratio between the class and the analytical result
ax2.semilogx(a,D/D_hypergeom(a,cosmo['LCDM']),lw=lw,label=r'$D_+/D_+^\mathrm{analytic}$')
ax2.semilogx(a,f/f_hypergeom(a,cosmo['LCDM']),lw=lw,label=r'$f/f^{\,\mathrm{analytic}}$')

# Draw a few vertical lines to orient ourselves
draw_vertical_redshift(cosmo['LCDM'], ax1, var='a',z=99,label='$z=99$')
draw_vertical_redshift(cosmo['LCDM'], ax1, var='a',z=49,label='$z=49$',ls='-')
draw_vertical_redshift(cosmo['LCDM'], ax2, var='a',z=99,label=None)
draw_vertical_redshift(cosmo['LCDM'], ax2, var='a',z=49,label=None,ls='-')

# Put some legend as well just for niceness
lgd1 = ax1.legend(fontsize=fs,ncol=1,loc='upper left',
           bbox_to_anchor=(1.02, 1.035))
lgd2 = ax2.legend(fontsize=fs,ncol=1,loc='upper left',
           bbox_to_anchor=(1.02, 0.83))

ax1.set_xlim([10**aminexp,1])
ax2.set_xlabel(r'$a$',fontsize=fs)
ax1.set_ylim([ymin,ymax])
ax2.set_ylim([0.9,1.099])

ax2.axhline(1,color='k')


fig.tight_layout()
fig.subplots_adjust(hspace=0.0)
fig.savefig('NewtonianGrowthFactor.pdf',bbox_extra_artists=(lgd1,lgd2), bbox_inches='tight')


# Possibly change the line properties
lw=2
fs=14
# Generate a new plot
fig, (ax1, ax2) = plt.subplots(2,1,figsize=(6,8),sharex=True,)
#                              gridspec_kw = {'height_ratios':[2, 1]})

# Now iterate through all cosmologies (with wa=0), not just LCDM
# for each one:
for M, csm in iter(cosmo.items()):
    # Get the w0/wa values for later labels
    if M!='LCDM':
        w0, wa = M.strip('()').split(',')
        # Let's ignore wa=/=0 for now
        if float(wa)!=0.0:
            continue

    # Then, for each case, get the class results
    bg = csm.get_background()
    a = 1./(bg['z']+1)
    H = bg['H [1/Mpc]']
    D = bg['gr.fac. D']
    f = bg['gr.fac. f']

    # Compare the class result to a bunch of approximations
    p=ax1.semilogx(a,D_linder2(a,csm)/a,lw=lw,ls='--',label=(r"$w_0={},~w_a={}$".format(w0,wa) if M!='LCDM' else M))
    # Repeat the color from the first line for all other lines
    color = p[0].get_color()
    ax1.semilogx(a,D/a,lw=lw,ls='-',color=color) # This is the class result, put as a solid line
    ax1.semilogx(a,D_integral2(a,csm)/a,lw=lw,ls=':',color=color)

    # For the other axes, plot the ratios of the class result over the given approximation
    ax2.semilogx(a,D/D_integral2(a,csm),lw=lw,ls='-',color=color)
    ax2.semilogx(a,D/D_hypergeom(a,csm),lw=lw,ls=':',color=color)
    ax2.semilogx(a,D/D_linder2(a,csm),lw=lw,ls='--',color=color)

#Just for getting labels for different line styles
ax1.plot([],[],color="grey",ls="--",label="Linder approximation")
ax1.plot([],[],color="grey",ls="-",label="CLASS")
ax1.plot([],[],color="grey",ls=":",label="Integral approximation")

ax2.plot([],[],color="grey",ls="-",label="Integral approximation")
ax2.plot([],[],color="grey",ls=":",label="Analytical approximation (only LCDM)")
ax2.plot([],[],color="grey",ls="--",label="Linder approximation")

ax1.set_xlim([1e-3,1])
ax2.set_xlabel(r'Scale factor $a$',fontsize=fs)
ax1.set_ylim([0,2])
ax2.set_ylim([0.9,1.3])
ax1.set_ylabel(r"$D(a)/a$")
ax2.set_ylabel(r"$D(a)/D^\mathrm{approx}(a)$")

lgd1 = ax1.legend(fontsize=fs,ncol=1,loc='lower left')
#           bbox_to_anchor=(1.0, 1.035))
lgd2 = ax2.legend(fontsize=fs,ncol=1,loc='upper right')

fig.tight_layout()
fig.subplots_adjust(hspace=0.0)
fig.savefig('Growthrate_w0.pdf')


# For those confused why the D_integral is performing so 'badly' (up to 10% error):
# The formula is not actually correct in dark energy scenarios with arbitrary equation of state w_DE(z)
# The true equation is the equation for D shown below, and we can see that numerically indeed it is very close to 0
# Transformed into Y = D/H, we can
csm = list(cosmo.values())[0]
bg = csm.get_background()
print(csm.pars)
a = 1./(bg['z']+1)
H = bg['H [1/Mpc]']
D = bg['gr.fac. D']

# check equation for D (with alpha=dlnH/dlna)
# D'' + D' (2+alpha) - 3/2 * Omega_m(a) * D === 0
Dp = CubicSpline(np.log(a),D).derivative()(np.log(a))
Dpp= CubicSpline(np.log(a),Dp).derivative()(np.log(a))
alpha = CubicSpline(np.log(a),np.log(H)).derivative()(np.log(a))
plt.semilogx(a,Dpp+Dp*(alpha+2)-1.5*(bg['(.)rho_cdm']+bg['(.)rho_b'])*D/H**2)

# Check equation for Y = D/H
# Y'' + Y' (3*alpha+2) + (2*alpha^2 + alpha' + 2*alpha - 3/2 * Omega_m(a)) * Y === 0
Y = D/H
Yp = CubicSpline(np.log(a),Y).derivative()(np.log(a))
Ypp = CubicSpline(np.log(a),Yp).derivative()(np.log(a))
alphaP = CubicSpline(np.log(a),alpha).derivative()(np.log(a))
plt.semilogx(a,Ypp+Yp*(3*alpha+2)+(2*alpha**2+alphaP+2*alpha-1.5*(bg['(.)rho_cdm']+bg['(.)rho_b'])/H**2)*Y)

# Finally, IFF 2*alpha^2 + alpha' + 2*alpha - 3/2 * Omega_m(a) == 0, then we would simply have
# Y'' + Y' (3*alpha+2) = 0, which would simplify for U = Y' to U' + U (3*alpha+2) = 0
# The solution to this is simply U = const * exp(-integral (3*alpha+2) dlna) = exp(-integral (3*dlnH/dlna+2)dlna) = exp(-3 ln H - 2 ln a) = 1/(a^2H^3)
# Now Y'=U => Y = integral U dlna = integral 1/(a^2 H^3) da/a = integral da/(a^3*H^3)
# So the usual integral expression holds ONLY IFF 2*alpha^2 + alpha' + 2*alpha - 3/2 * Omega_m(a) == 0
# But, it can be shown that for general w_DE that only holds if w_DE(a) = -1 or w_DE(a)=-1/3
# For this, you can use that iff w_DE is constant that
# alpha = -3/2 * (1+w_DE * (1-Omega_m(a)))
# alpha' = 9/2 * w_DE^2 * Omega_m(a) * (1-Omega_m(a))
# Plugging everything in, with w_DE = -1 then leads to the desired canclation
# However, it is no wonder that in w_DE = -0.7, for example, the integral approximation is not great, since there is no reason for it to hold!
