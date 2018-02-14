import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/net/home/lxtsfs1/tpe/archidiacono/codes/class_devel/output/lcdm_cl_lensed.dat', '/net/home/lxtsfs1/tpe/archidiacono/codes/class_devel/output/nuint0_cl_lensed.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
roots = ['lcdm_cl_lensed', 'nuint0_cl_lensed']

fig, ax = plt.subplots()
y_axis = [u'TT']
tex_names = ['TT']
x_axis = 'l'
ax.set_xlabel('$\ell$', fontsize=16)
plt.show()
fig.savefig('nuint.eps')