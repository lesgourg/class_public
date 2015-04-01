import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/Users/poulin/Documents/Doc_Labo/ProgrammeDarkAges/class_public-2.4.2/DM_Annihilation/fEM_muons_1GeV_smallhalos.dat', '/Users/poulin/Documents/Doc_Labo/ProgrammeDarkAges/class_public-2.4.2/DM_Annihilation/fEM_electron_1GeV_smallhalos.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
roots = ['fEM_muons_1GeV_smallhalos', 'fEM_electron_1GeV_smallhalos']

fig, ax = plt.subplots()

index, curve = 0, data[0]
y_axis = [u'fEM']
tex_names = ['fEM']
x_axis = 'z'
ylim = []
xlim = [10.0, 10000.0]
ax.semilogx(curve[:, 0], curve[:, 1])

index, curve = 1, data[1]
y_axis = [u'fEM']
tex_names = ['fEM']
x_axis = 'z'
ylim = []
xlim = [10.0, 10000.0]
ax.semilogx(curve[:, 0], curve[:, 1])

ax.legend([root+': '+elem for (root, elem) in
    itertools.product(roots, y_axis)], loc='best')

ax.set_xlabel('z', fontsize=16)
ax.set_xlim(xlim)
ax.set_ylim()
plt.show()