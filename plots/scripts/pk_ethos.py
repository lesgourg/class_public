import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/home/archidiacono/class_public-2.4.3_ethos/output/lcdm_pk.dat', '/home/archidiacono/class_public-2.4.3_ethos/output/boehm5_pk.dat','/home/archidiacono/class_public-2.4.3_ethos/output/boehm6_pk.dat', '/home/archidiacono/class_public-2.4.3_ethos/output/boehm7_pk.dat', '/home/archidiacono/class_public-2.4.3_ethos/output/boehm8_pk.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
#roots = ['nudm3_prec_000_cl_lensed', 'nudm3_prec_001_cl_lensed']

fig, ax = plt.subplots()

index, curve = 0, data[0]
y_axis = [u'Pk']
tex_names = ['Pk']
x_axis = 'k'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 1]), linestyle='-', color='k')

index, curve = 1, data[1]
y_axis = [u'Pk']
tex_names = ['Pk']
x_axis = 'k'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 1]),linestyle=':',color='b')

index, curve = 1, data[2]
y_axis = [u'Pk']
tex_names = ['Pk']
x_axis = 'k'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 1]),linestyle=':',color='g')

index, curve = 1, data[3]
y_axis = [u'Pk']
tex_names = ['Pk']
x_axis = 'k'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 1]),linestyle='--',color='r')

index, curve = 1, data[4]
y_axis = [u'Pk']
tex_names = ['Pk']
x_axis = 'k'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 1]),linestyle='-.',color='r')

ax.set_xlim([0.3,200])
ax.set_ylim([1.e-10,1.e4])

ax.legend(['$u=0$','$u=10^{-5}$', '$u=10^{-6}$', '$u=10^{-7}$', '$u=10^{-8}$'], loc='best')

ax.set_xlabel('$k [h/Mpc]$', fontsize=16)
ax.set_ylabel('$P(k) \, [(Mpc/h)^3]$', fontsize=16)
plt.show()
plt.savefig('pk_boehm.pdf')
