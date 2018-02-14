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

files=['lcdm_z2_pk_nl.dat','pm_ncdm_adjOc_z2_pk_nl.dat','pcb_ncdm_adjOc_z2_pk_nl.dat','pm_ncdm_adjOc_z2_tk.dat', 'Bird_pcb_ncdm_adjOc_z2_pk_nl.dat'] #lcdm00,mcdm01,mcdm(no neutrino transfer)02

data = []
for data_file in files:
    data.append(np.loadtxt(data_file))

fig, ax = plt.subplots()

lcdm = data[0]
m = data[1]
cb = data[2]
transf=data[3]
tk_c=transf[:,3]
tk_b=transf[:,2]
tk_nu=transf[:,5]
Bird=data[4]

interpolation = InterpolatedUnivariateSpline(m[:,0],m[:,1])
pk_m = interpolation(lcdm[:,0])
ax.plot(lcdm[:,0], pk_m[:]/lcdm[:,1]-1., linestyle='-', color='b')

interpolation = InterpolatedUnivariateSpline(cb[:,0],cb[:,1])
pk_cb = interpolation(lcdm[:,0])
ax.plot(lcdm[:,0], pk_cb[:]/lcdm[:,1]-1., linestyle='--', color='b')

interpolation = InterpolatedUnivariateSpline(transf[:,0],tk_c[:])
tk_c = interpolation(lcdm[:,0])
interpolation = InterpolatedUnivariateSpline(transf[:,0],tk_b[:])
tk_b = interpolation(lcdm[:,0])
interpolation = InterpolatedUnivariateSpline(transf[:,0],tk_nu[:])
tk_nu = interpolation(lcdm[:,0])
tk_bc=(oc*tk_c+ob*tk_b)/(oc+ob)
tk_m=(oc*tk_c+ob*tk_b+Mnu/93.14*tk_nu)/(oc+ob+Mnu/93.14)
pk_cb_from_transf=(tk_bc/tk_m)**2*pk_m
ax.plot(lcdm[:,0], pk_cb_from_transf[:]/lcdm[:,1]-1., linestyle=':', color='g')

interpolation = InterpolatedUnivariateSpline(Bird[:,0],Bird[:,1])
pk_cb_Bird = interpolation(lcdm[:,0])
ax.plot(lcdm[:,0], pk_cb_Bird[:]/lcdm[:,1]-1., linestyle='--', color='r')

ax.set_xscale('log')
#ax.set_xlim([0.002,0.6])
#ax.set_ylim([-0.1,0.1])

plt.axhline(-8.*Mnu/93.14/(oc+ob+Mnu/93.14),color='k', linestyle=':', linewidth=0.1)
ax.legend([r'full (Bird&Viel)', r'CDM+baryons (Takahashi)', r'full (Bird&Viel) x T(k) scaling', r'CDM+baryons (Bird&Viel)'],loc='best',fontsize=12,fancybox=True, framealpha=0.5)
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

