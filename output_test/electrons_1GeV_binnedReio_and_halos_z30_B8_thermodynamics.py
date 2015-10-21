import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/Users/poulin/Documents/Doc_Labo/ProgrammeDarkAges/class_public-2.4.2/output_test/electrons_1GeV_binnedReio_and_halos_z30_B8_thermodynamics.dat', '/Users/poulin/Documents/Doc_Labo/ProgrammeDarkAges/class_public-2.4.2/output_test/electrons_1GeV_binnedReio_thermodynamics.dat', '/Users/poulin/Documents/Doc_Labo/ProgrammeDarkAges/class_public-2.4.2/output_cl/Muons_1GeV_Halo_zh30_B8_0079_thermodynamics.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
roots = ['electrons_1GeV_binnedReio_and_halos_z30_B8_thermodynamics', 'electrons_1GeV_binnedReio_thermodynamics', 'Muons_1GeV_Halo_zh30_B8_0079_thermodynamics']

fig, ax = plt.subplots()

index, curve = 0, data[0]
y_axis = [u'x_e']
tex_names = ['x_e']
x_axis = 'z'
ylim = []
xlim = [5.0, 12.0]
ax.plot(curve[:, 0], curve[:, 2])

index, curve = 1, data[1]
y_axis = [u'x_e']
tex_names = ['x_e']
x_axis = 'z'
ylim = []
xlim = [5.0, 12.0]
ax.plot(curve[:, 0], curve[:, 2])

index, curve = 2, data[2]
y_axis = [u'x_e']
tex_names = ['x_e']
x_axis = 'z'
ylim = []
xlim = [5.0, 12.0]
ax.plot(curve[:, 0], curve[:, 2])

ax.legend([root+': '+elem for (root, elem) in
    itertools.product(roots, y_axis)], loc='best')

ax.set_xlabel('z', fontsize=16)
ax.set_xlim(xlim)
ax.set_ylim()
plt.show()