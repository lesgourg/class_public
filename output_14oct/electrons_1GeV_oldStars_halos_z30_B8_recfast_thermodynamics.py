import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/Users/poulin/Documents/Doc_Labo/ProgrammeDarkAges/class_public-2.4.2/output_14oct/electrons_1GeV_oldStars_halos_z30_B8_recfast_thermodynamics.dat', '/Users/poulin/Documents/Doc_Labo/ProgrammeDarkAges/class_public-2.4.2/output_14oct/electrons_1GeV_oldStars_halos_z30_B8_thermodynamics.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
roots = ['electrons_1GeV_oldStars_halos_z30_B8_recfast_thermodynamics', 'electrons_1GeV_oldStars_halos_z30_B8_thermodynamics']

fig, ax = plt.subplots()

index, curve = 0, data[0]
y_axis = [u'Tb[K]']
tex_names = ['Tb [K]']
x_axis = 'z'
ylim = []
xlim = [1.0, 300.0]
ax.loglog(curve[:, 0], abs(curve[:, 6]))

index, curve = 1, data[1]
y_axis = [u'Tb[K]']
tex_names = ['Tb [K]']
x_axis = 'z'
ylim = []
xlim = [1.0, 300.0]
ax.loglog(curve[:, 0], abs(curve[:, 6]))

ax.legend([root+': '+elem for (root, elem) in
    itertools.product(roots, y_axis)], loc='best')

ax.set_xlabel('z', fontsize=16)
ax.set_xlim(xlim)
ax.set_ylim()
plt.show()