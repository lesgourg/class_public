import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/net/home/lxtsfs1/tpe/archidiacono/codes/class/output/pk05_cl_lensed.dat', '/net/home/lxtsfs1/tpe/archidiacono/codes/class/output/pk06_cl_lensed.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
roots = ['pk05_cl_lensed', 'pk06_cl_lensed']

fig, ax = plt.subplots()
y_axis = [u'phiphi']
tex_names = ['phiphi']
x_axis = 'l'
ax.set_xlabel('$\ell$', fontsize=16)
plt.show()
fig.savefig('tt.pdf')