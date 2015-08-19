import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/Users/poulin/Documents/Doc_Labo/ProgrammeDarkAges/class_public-2.4.2/output_test/Reio_DM_electrons_1GeV_B12_zh30_thermodynamics.dat', '/Users/poulin/Documents/Doc_Labo/ProgrammeDarkAges/class_public-2.4.2/output_test/Reio_Stars_thermodynamics.dat', '/Users/poulin/Documents/Doc_Labo/ProgrammeDarkAges/class_public-2.4.2/output_test/Reio_None_thermodynamics.dat', '/Users/poulin/Documents/Doc_Labo/ProgrammeDarkAges/class_public-2.4.2/output_test/Reio_DM_No_halos_thermodynamics.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
roots = ['Reio_DM_electrons_1GeV_B12_zh30_thermodynamics', 'Reio_Stars_thermodynamics', 'Reio_None_thermodynamics', 'Reio_DM_No_halos_thermodynamics']

fig, ax = plt.subplots()

index, curve = 0, data[0]
y_axis = [u'g[Mpc^-1]']
tex_names = ['g [Mpc^-1]']
x_axis = 'z'
ylim = [1e-11, 1.0]
xlim = [1.0, 1300.0]
ax.loglog(curve[:, 0], abs(curve[:, 5]))

index, curve = 1, data[1]
y_axis = [u'g[Mpc^-1]']
tex_names = ['g [Mpc^-1]']
x_axis = 'z'
ylim = [1e-11, 1.0]
xlim = [1.0, 1300.0]
ax.loglog(curve[:, 0], abs(curve[:, 5]))

index, curve = 2, data[2]
y_axis = [u'g[Mpc^-1]']
tex_names = ['g [Mpc^-1]']
x_axis = 'z'
ylim = [1e-11, 1.0]
xlim = [1.0, 1300.0]
ax.loglog(curve[:, 0], abs(curve[:, 5]))

index, curve = 3, data[3]
y_axis = [u'g[Mpc^-1]']
tex_names = ['g [Mpc^-1]']
x_axis = 'z'
ylim = [1e-11, 1.0]
xlim = [1.0, 1300.0]
ax.loglog(curve[:, 0], abs(curve[:, 5]))

ax.legend([root+': '+elem for (root, elem) in
    itertools.product(roots, y_axis)], loc='best')

ax.set_xlabel('z', fontsize=16)
ax.set_xlim(xlim)
ax.set_ylim()
ax.set_ylim(ylim)
plt.show()