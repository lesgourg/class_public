import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/net/home/lxtsfs1/tpe/archidiacono/codes/class/output/lcdm_new_z1_pk.dat', '/net/home/lxtsfs1/tpe/archidiacono/codes/class/output/lcdm_syn_z1_pk.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
roots = ['lcdm_new_z1_pk', 'lcdm_syn_z1_pk']

fig, ax = plt.subplots()
y_axis = [u'P(Mpc/h)^3']
tex_names = ['P (Mpc/h)^3']
x_axis = 'k (h/Mpc)'
ax.set_xlabel('k (h/Mpc)', fontsize=16)
plt.show()
fig.savefig('lcdm_syn_vs_new.pdf')