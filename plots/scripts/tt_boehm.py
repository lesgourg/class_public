import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/home/archidiacono/ethos_class/output/lcdm_cl_lensed.dat', '/home/archidiacono/ethos_class/output/boehm_T2_u1_cl_lensed.dat','/home/archidiacono/ethos_class/output/boehm_T2_u2_cl_lensed.dat', '/home/archidiacono/ethos_class/output/boehm_T2_u3_cl_lensed.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
#roots = ['nudm3_prec_000_cl_lensed', 'nudm3_prec_001_cl_lensed']

fig, ax = plt.subplots()

index, curve = 0, data[0]
y_axis = [u'TT']
tex_names = ['TT']
x_axis = 'l'
ylim = []
xlim = []
ax.plot(curve[:, 0], abs(curve[:, 1])*7.42835025e12, linestyle='-', color='k')

index, curve = 1, data[1]
y_axis = [u'TT']
tex_names = ['TT']
x_axis = 'l'
ylim = []
xlim = []
ax.plot(curve[:, 0], abs(curve[:, 1])*7.42835025e12,linestyle='-.',color='r')

index, curve = 1, data[2]
y_axis = [u'TT']
tex_names = ['TT']
x_axis = 'l'
ylim = []
xlim = []
ax.plot(curve[:, 0], abs(curve[:, 1])*7.42835025e12,linestyle=':',color='g')

index, curve = 1, data[3]
y_axis = [u'TT']
tex_names = ['TT']
x_axis = 'l'
ylim = []
xlim = []
ax.plot(curve[:, 0], abs(curve[:, 1])*7.42835025e12,linestyle='--',color='r')

ax.legend(['u=0','u=0.1', 'u=0.01', 'u=0.001'], loc='best')

ax.set_xlabel('$\ell$', fontsize=16)
ax.set_ylabel('$\ell(\ell+1)C_{\ell}^{TT}/(2\pi) \, [\mu K^2]$', fontsize=16)
plt.show()
plt.savefig('tt_boehm.pdf')
