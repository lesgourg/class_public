import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/Users/poulin/Documents/Doc_Labo/ProgrammeDarkAges/class_public-2.4.2/output_14oct/electrons_1GeV_ModelStarsReio_and_halos_z30_B8_thermodynamics.dat', '/Users/poulin/Documents/Doc_Labo/ProgrammeDarkAges/class_public-2.4.2/output_14oct/ModelStarsReio_AvecHelium_thermodynamics.dat', '/Users/poulin/Documents/Doc_Labo/ProgrammeDarkAges/class_public-2.4.2/output_23oct/Standard_planck_zreio65_thermodynamics.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
roots = ['electrons_1GeV_ModelStarsReio_and_halos_z30_B8_thermodynamics', 'ModelStarsReio_AvecHelium_thermodynamics', 'Standard_planck_zreio65_thermodynamics']

fig, ax = plt.subplots()

index, curve = 0, data[0]
y_axis = [u'g[Mpc^-1]']
tex_names = ['g [Mpc^-1]']
x_axis = 'z'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 5]))

index, curve = 1, data[1]
y_axis = [u'g[Mpc^-1]']
tex_names = ['g [Mpc^-1]']
x_axis = 'z'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 5]))

index, curve = 2, data[2]
y_axis = [u'g[Mpc^-1]']
tex_names = ['g [Mpc^-1]']
x_axis = 'z'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 5]))

ax.legend([root+': '+elem for (root, elem) in
    itertools.product(roots, y_axis)], loc='best')

ax.set_xlabel('z', fontsize=16)
plt.show()