import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/Users/poulin/Documents/Doc_Labo/ProgrammeDarkAges/class_public-2.4.2/output_new/test_thermodynamics.dat', '/Users/poulin/Documents/Doc_Labo/ProgrammeDarkAges/class_public-2.4.2/output_new/NoReio_thermodynamics.dat', '/Users/poulin/Documents/Doc_Labo/ProgrammeDarkAges/class_public-2.4.2/output_new/electrons_1GeV_zh30_B8_taureio79_thermodynamics.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
roots = ['test_thermodynamics', 'NoReio_thermodynamics', 'electrons_1GeV_zh30_B8_taureio79_thermodynamics']

fig, ax = plt.subplots()

index, curve = 0, data[0]
y_axis = [u'g[Mpc^-1]']
tex_names = ['g [Mpc^-1]']
x_axis = 'conf. time [Mpc]'
ylim = [1e-21, 1.0]
xlim = []
ax.loglog(curve[:, 1], abs(curve[:, 5]))

index, curve = 1, data[1]
y_axis = [u'g[Mpc^-1]']
tex_names = ['g [Mpc^-1]']
x_axis = 'conf. time [Mpc]'
ylim = [1e-21, 1.0]
xlim = []
ax.loglog(curve[:, 1], abs(curve[:, 5]))

index, curve = 2, data[2]
y_axis = [u'g[Mpc^-1]']
tex_names = ['g [Mpc^-1]']
x_axis = 'conf. time [Mpc]'
ylim = [1e-21, 1.0]
xlim = []
ax.loglog(curve[:, 1], abs(curve[:, 5]))

ax.legend([root+': '+elem for (root, elem) in
    itertools.product(roots, y_axis)], loc='best')

ax.set_xlabel('conf. time [Mpc]', fontsize=16)
ax.set_ylim(ylim)
plt.show()