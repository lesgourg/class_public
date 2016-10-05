import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/home/archidiacono/ethos_class/output/lcdm_cl_lensed.dat', '/home/archidiacono/ethos_class/output/ethos_fig1_n1_cl_lensed.dat','/home/archidiacono/ethos_class/outcamb/ethos_fig1_n1_lensedCls.dat',
'/home/archidiacono/ethos_class/output/ethos_fig1_n2_cl_lensed.dat','/home/archidiacono/ethos_class/outcamb/ethos_fig1_n2_lensedCls.dat',
'/home/archidiacono/ethos_class/output/ethos_fig1_n3_cl_lensed.dat','/home/archidiacono/ethos_class/outcamb/ethos_fig1_n3_lensedCls.dat',
'/home/archidiacono/ethos_class/output/ethos_fig1_n4_cl_lensed.dat','/home/archidiacono/ethos_class/outcamb/ethos_fig1_n4_lensedCls.dat']


data = []
for data_file in files:
    data.append(np.loadtxt(data_file))

fig, ax = plt.subplots()

index, curve = 0, data[0]
y_axis = [u'TT']
tex_names = ['TT']
x_axis = 'l'
ylim = []
xlim = []
ax.plot(curve[:, 0], abs(curve[:, 1])*7.42835025e12, linestyle=':', color='k')

index, curve = 1, data[1]
y_axis = [u'TT']
tex_names = ['TT']
x_axis = 'l'
ylim = []
xlim = []
ax.plot(curve[:, 0], abs(curve[:, 1])*7.42835025e12,linestyle='--',color='r')

index, curve = 2, data[2]
y_axis = [u'TT']
tex_names = ['TT']
x_axis = 'l'
ylim = []
xlim = []
ax.plot(curve[:, 0], abs(curve[:, 1]),linestyle='-',color='r')

index, curve = 3, data[3]
y_axis = [u'TT']
tex_names = ['TT']
x_axis = 'l'
ylim = []
xlim = []
ax.plot(curve[:, 0], abs(curve[:, 1])*7.42835025e12,linestyle='--',color='b')

index, curve = 4, data[4]
y_axis = [u'TT']
tex_names = ['TT']
x_axis = 'l'
ylim = []
xlim = []
ax.plot(curve[:, 0], abs(curve[:, 1]), linestyle='-', color='b')

index, curve = 5, data[5]
y_axis = [u'TT']
tex_names = ['TT']
x_axis = 'l'
ylim = []
xlim = []
ax.plot(curve[:, 0], abs(curve[:, 1])*7.42835025e12,linestyle='--',color='g')

index, curve = 6, data[6]
y_axis = [u'TT']
tex_names = ['TT']
x_axis = 'l'
ylim = []
xlim = []
ax.plot(curve[:, 0], abs(curve[:, 1]),linestyle='-',color='g')

index, curve = 7, data[7]
y_axis = [u'TT']
tex_names = ['TT']
x_axis = 'l'
ylim = []
xlim = []
ax.plot(curve[:, 0], abs(curve[:, 1])*7.42835025e12,linestyle='--',color='k')

index, curve = 8, data[8]
y_axis = [u'TT']
tex_names = ['TT']
x_axis = 'l'
ylim = []
xlim = []
ax.plot(curve[:, 0], abs(curve[:, 1]),linestyle='-',color='k')

ax.legend(['class n=1','camb n=1','class n=2','camb n=2','class n=3','camb n=3','class n=4','camb n=4'], loc='best')

ax.set_xlabel('$\ell$', fontsize=16)
ax.set_ylabel('$\ell(\ell+1)C_{\ell}^{TT}/(2\pi) \, [\mu K^2]$', fontsize=16)
plt.show()
plt.savefig('tt_ethos.pdf')
