import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import itertools
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.interpolate import interp1d


myz=1.0
oc=0.12038
ob=0.022032
Mnu=3.*0.02

files=['pk00_z2_pk.dat','pk05_z2_pk.dat','pk06_z2_pk.dat','pk05_z2_tk.dat'] #lcdm00,mcdm01,mcdm(no neutrino transfer)02

data = []
for data_file in files:
    data.append(np.loadtxt(data_file))

fig, ax = plt.subplots()

lcdm = data[0]
ncdm = data[1]
ncdm_bc = data[2]
ncdm_transf=data[3]
tk_c=ncdm_transf[:,3]
tk_b=ncdm_transf[:,2]
tk_nu=ncdm_transf[:,5]

interpolation = InterpolatedUnivariateSpline(ncdm[:,0],ncdm[:,1])
pk_ncdm = interpolation(lcdm[:,0])
ax.plot(lcdm[:,0], pk_ncdm[:]/lcdm[:,1]-1., linestyle='-', color='b')

interpolation = InterpolatedUnivariateSpline(ncdm_bc[:,0],ncdm_bc[:,1])
interpolated = interpolation(lcdm[:,0])
ax.plot(lcdm[:,0], interpolated[:]/lcdm[:,1]-1., linestyle='--', color='b')

interpolation = InterpolatedUnivariateSpline(ncdm_transf[:,0],tk_c[:])
tk_c = interpolation(lcdm[:,0])
interpolation = InterpolatedUnivariateSpline(ncdm_transf[:,0],tk_b[:])
tk_b = interpolation(lcdm[:,0])
interpolation = InterpolatedUnivariateSpline(ncdm_transf[:,0],tk_nu[:])
tk_nu = interpolation(lcdm[:,0])
tk_bc=(oc*tk_c+ob*tk_b)/(oc+ob)
tk_m=(oc*tk_c+ob*tk_b+Mnu/93.14*tk_nu)/(oc+ob+Mnu/93.14)
pk_ncdm_cb=(tk_bc/tk_m)**2*pk_ncdm
ax.plot(lcdm[:,0], pk_ncdm_cb[:]/lcdm[:,1]-1., linestyle=':', color='g')

ax.set_xscale('log')
#ax.set_xlim([0.002,0.6])
#ax.set_ylim([-0.1,0.1])

plt.axhline(-8.*Mnu/93.14/(oc+ob+Mnu/93.14),color='k', linestyle=':', linewidth=0.1)
ax.legend([r'full',r'CDM+baryons', r'CDM+baryons, T(k) scaling'],loc='best',fontsize=12,fancybox=True, framealpha=0.5)
#ax.legend(handles=[varm,adjtheta,adjevery,desi],loc='best')

textstr = '$z=%.1f$'%(myz)
ax.text(0.05, 0.97, textstr, transform=ax.transAxes, fontsize=12,
        verticalalignment='top')#, bbox=props)

#ax.tick_params(axis='x', labelsize=14)
#ax.tick_params(axis='y', labelsize=14)

pos1 = ax.get_position() # get the original position 
pos2 = [pos1.x0 + 0.01, pos1.y0,  pos1.width, pos1.height] 
ax.set_position(pos2) # set a new position

ax.set_ylabel(r'$P^{M_\nu}(k)/P^{M_\nu=0}(k)-1$', fontsize=12,labelpad=0)
ax.set_xlabel('$k\,[h/\mathrm{Mpc}]$', fontsize=12,labelpad=0)
plt.show()
plt.savefig('diff_pk_plot.pdf')

