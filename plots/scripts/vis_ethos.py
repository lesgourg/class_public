import matplotlib.pyplot as plt
import numpy as np
import itertools
from scipy.interpolate import InterpolatedUnivariateSpline

files = ['/home/archidiacono/ethos_class/output/ethos_fig1_n2_thermodynamics.dat','/home/archidiacono/ethos_class/output/ethos_fig1_n2_N3046_thermodynamics.dat','/home/archidiacono/ethos_class/output/ethos_fig1_n4_thermodynamics.dat','/home/archidiacono/ethos_class/output/ethos_fig1_n4_N3046_thermodynamics.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
#roots = ['nudm3_prec_000_cl_lensed', 'nudm3_prec_001_cl_lensed']

fig, ax = plt.subplots()

ax.set_xscale('log')

index, thermo = 0, data[0]
ax.plot(thermo[:, 1], thermo[:,13], linestyle='-', color='b')

index, thermo = 1, data[1]
ax.plot(thermo[:, 1], thermo[:,13], linestyle='--', color='b')

index, thermo = 2, data[2]
ax.plot(thermo[:, 1], thermo[:,13], linestyle='-', color='k')

index, thermo = 3, data[3]
ax.plot(thermo[:, 1], thermo[:,13], linestyle='--', color='k')

ax.set_xlim([5.e-3,1.e1])

ax.legend(['Nur=3.046 Nint=0.25 n=2', 'Nur=0 Nint=3.046 n=2','Nur=3.046 Nint=0.25 n=4','Nur=0 Nint=3.046 n=4'], loc='best')

ax.set_xlabel('$tau$ Mpc', fontsize=16)
ax.set_ylabel('$dmu \mathrm{exp}(-mu) \mathrm{Mpc}^{-1}$', fontsize=16)
plt.show()
plt.savefig('vis_ethos_N3046.pdf')
