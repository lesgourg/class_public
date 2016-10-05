import matplotlib.pyplot as plt
import numpy as np
import itertools
from scipy.interpolate import InterpolatedUnivariateSpline

files = ['/home/archidiacono/ethos_class/output/lcdm_cl_lensed.dat', '/home/archidiacono/ethos_class/output/test_N3_T2_a1e5_s1_cl_lensed.dat','/home/archidiacono/ethos_class/output/test_N3_T2_a1e5_s0_cl_lensed.dat']


data = []
for data_file in files:
    data.append(np.loadtxt(data_file))

fig, ax = plt.subplots()

index, curve0 = 0, data[0]
#ax.plot(curve[:, 0], abs(curve[:, 1])*7.42835025e12, linestyle=':', color='k')

index, curve = 1, data[1]
ax.plot(curve0[:, 0], curve[:, 1]/curve0[:,1]-1.,linestyle='-',color='k')

index, curve = 1, data[2]
ax.plot(curve0[:, 0], curve[:, 1]/curve0[:,1]-1.,linestyle='--',color='k')

ax.legend(['Nur=0 Nint=3.046 n=2 s=1','Nur=0 Nint=3.046 n=2 s=0'], loc='best')

ax.set_xlabel('$\ell$', fontsize=16)
ax.set_ylabel('$C_{\ell}/C_{\ell}^{LCDM}-1$', fontsize=16)
plt.show()
plt.savefig('tt_err_test_N3_T2_a1e5.pdf')
