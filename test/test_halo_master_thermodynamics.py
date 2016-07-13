import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/Users/poulin/Documents/Labo/ProgrammeDarkAges/class_public-2.4.2/test/test_halo_master_thermodynamics.dat', '/Users/poulin/Documents/Labo/ProgrammeDarkAges/class_public-2.4.2/test/test_halo_master_beyond_thermodynamics.dat', '/Users/poulin/Documents/Labo/ProgrammeDarkAges/class_public-2.4.2/test/test_halo_perso_thermodynamics.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
roots = ['test_halo_master_thermodynamics', 'test_halo_master_beyond_thermodynamics', 'test_halo_perso_thermodynamics']

fig, ax = plt.subplots()

index, curve = 0, data[0]
y_axis = [u'x_e']
tex_names = ['x_e']
x_axis = 'z'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 2]))

index, curve = 1, data[1]
y_axis = [u'x_e']
tex_names = ['x_e']
x_axis = 'z'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 2]))

index, curve = 2, data[2]
y_axis = [u'x_e']
tex_names = ['x_e']
x_axis = 'z'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 2]))

ax.legend([root+': '+elem for (root, elem) in
    itertools.product(roots, y_axis)], loc='best')

ax.set_xlabel('z', fontsize=16)
plt.show()