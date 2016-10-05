import matplotlib.pyplot as plt
import numpy as np
import itertools
from scipy.interpolate import InterpolatedUnivariateSpline

plt.rc('text', usetex=True)

files = ['/home/archidiacono/ethos_class/output/lcdm_pk.dat', '/home/archidiacono/ethos_class/output/ethos_fig1_n1_pk.dat','/home/archidiacono/ethos_class/outcamb/ethos_fig1_n1_matterpower.dat',
'/home/archidiacono/ethos_class/output/ethos_fig1_n2_pk.dat','/home/archidiacono/ethos_class/outcamb/ethos_fig1_n2_matterpower.dat',
'/home/archidiacono/ethos_class/output/ethos_fig1_n3_pk.dat','/home/archidiacono/ethos_class/outcamb/ethos_fig1_n3_matterpower.dat',
'/home/archidiacono/ethos_class/output/ethos_fig1_n4_pk.dat','/home/archidiacono/ethos_class/outcamb/ethos_fig1_n4_matterpower.dat']

data = []
for data_file in files:
    data.append(np.loadtxt(data_file))

fig, ax = plt.subplots()

ax.set_xscale('log')
ax.set_xlim([1,200])
ax.set_xlabel(r'$k \,\,\, [h/\mathrm{Mpc}]$', fontsize=16)

ax.set_yscale('log')
ax.set_ylim([1.e-12,10])
ax.set_ylabel(r'$P(k)~/~P(k) \, [\Lambda \mathrm{CDM}]$', fontsize=16)

ref = data[0]
#plt.bar(planck[:,1],planck[:,4]/planck[:,3],width=planck[:,2]-planck[:,1],bottom=-planck[:,4]/planck[:,3]/2.,color='lavenderblush',edgecolor='lightgray')

model = data[1]
interpolation = InterpolatedUnivariateSpline(model[:,0],model[:,1])
interpolated = interpolation(ref[:,0])
ax.plot(ref[:,0], interpolated[:]/ref[:,1], color='r', linestyle='--')

model = data[2]
interpolation = InterpolatedUnivariateSpline(model[:,0],model[:,1])
interpolated = interpolation(ref[:,0])
ax.plot(ref[:,0], interpolated[:]/ref[:,1], color='r')

model = data[3]
interpolation = InterpolatedUnivariateSpline(model[:,0],model[:,1])
interpolated = interpolation(ref[:,0])
ax.plot(ref[:,0], interpolated[:]/ref[:,1], color='b', linestyle='--')

model = data[4]
interpolation = InterpolatedUnivariateSpline(model[:,0],model[:,1])
interpolated = interpolation(ref[:,0])
ax.plot(ref[:,0], interpolated[:]/ref[:,1], color='b')

model = data[5]
interpolation = InterpolatedUnivariateSpline(model[:,0],model[:,1])
interpolated = interpolation(ref[:,0])
ax.plot(ref[:,0], interpolated[:]/ref[:,1], color='g', linestyle='--')

model = data[6]
interpolation = InterpolatedUnivariateSpline(model[:,0],model[:,1])
interpolated = interpolation(ref[:,0])
ax.plot(ref[:,0], interpolated[:]/ref[:,1], color='g')

model = data[7]
interpolation = InterpolatedUnivariateSpline(model[:,0],model[:,1])
interpolated = interpolation(ref[:,0])
ax.plot(ref[:,0], interpolated[:]/ref[:,1], color='k', linestyle='--')

model = data[8]
interpolation = InterpolatedUnivariateSpline(model[:,0],model[:,1])
interpolated = interpolation(ref[:,0])
ax.plot(ref[:,0], interpolated[:]/ref[:,1], color='k')

ax.legend(['class n=1','camb n=1','class n=2','camb n=2','class n=3','camb n=3','class n=4','camb n=4'], loc='best')

#plt.text(1.1e-4,0.05,'$N_{ur}=0$, $N_{ncdm}=1$, $m_{ncdm}=0.06$eV, $num-deg=3$, $N_{dg}=0.21$', ha='left', va='center')

#plt.axhline(0,color='k', linestyle='--')

#plt.show()
plt.savefig('pk_ratio_ethos_prova.pdf')
