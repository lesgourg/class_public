import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/Users/poulin/Documents/Doc_Labo/ProgrammeDarkAges/class_public-2.4.2/output_14oct/ModelStarsReio_AvecHelium_cl.dat', '/Users/poulin/Documents/Doc_Labo/ProgrammeDarkAges/class_public-2.4.2/output_14oct/electrons_1GeV_ModelStarsReio_DM_noHalos_cl.dat', '/Users/poulin/Documents/Doc_Labo/ProgrammeDarkAges/class_public-2.4.2/output_14oct/electrons_1GeV_ModelStarsReio_and_halos_z30_B8_cl.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
roots = ['ModelStarsReio_AvecHelium_cl', 'electrons_1GeV_ModelStarsReio_DM_noHalos_cl', 'electrons_1GeV_ModelStarsReio_and_halos_z30_B8_cl']

fig, ax = plt.subplots()

index, curve = 0, data[0]
y_axis = [u'TT', u'TE']
tex_names = ['TT', 'TE']
x_axis = 'l'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 1]))
ax.loglog(curve[:, 0], abs(curve[:, 3]))

index, curve = 1, data[1]
y_axis = [u'TT', u'TE']
tex_names = ['TT', 'TE']
x_axis = 'l'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 1]))
ax.loglog(curve[:, 0], abs(curve[:, 3]))

index, curve = 2, data[2]
y_axis = [u'TT', u'TE']
tex_names = ['TT', 'TE']
x_axis = 'l'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 1]))
ax.loglog(curve[:, 0], abs(curve[:, 3]))

ax.legend([root+': '+elem for (root, elem) in
    itertools.product(roots, y_axis)], loc='best')

ax.set_xlabel('$\ell$', fontsize=16)
plt.show()