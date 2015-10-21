import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/Users/poulin/Documents/Doc_Labo/ProgrammeDarkAges/class_public-2.4.2/output_14oct/Standard_Reio_0079_cl.dat', '/Users/poulin/Documents/Doc_Labo/ProgrammeDarkAges/class_public-2.4.2/output_14oct/ModelStarsReio_AvecHelium_cl.dat', '/Users/poulin/Documents/Doc_Labo/ProgrammeDarkAges/class_public-2.4.2/output_14oct/ModelStarsReio_SansHelium_cl.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
roots = ['Standard_Reio_0079_cl', 'ModelStarsReio_AvecHelium_cl', 'ModelStarsReio_SansHelium_cl']

fig, ax = plt.subplots()

index, curve = 0, data[0]
y_axis = [u'EE']
tex_names = ['EE']
x_axis = 'l'
ylim = []
xlim = [0.0, 2000.0]
ax.loglog(curve[:, 0], abs(curve[:, 2]))

index, curve = 1, data[1]
y_axis = [u'EE']
tex_names = ['EE']
x_axis = 'l'
ylim = []
xlim = [0.0, 2000.0]
ax.loglog(curve[:, 0], abs(curve[:, 2]))

index, curve = 2, data[2]
y_axis = [u'EE']
tex_names = ['EE']
x_axis = 'l'
ylim = []
xlim = [0.0, 2000.0]
ax.loglog(curve[:, 0], abs(curve[:, 2]))

ax.legend([root+': '+elem for (root, elem) in
    itertools.product(roots, y_axis)], loc='best')

ax.set_xlabel('$\ell$', fontsize=16)
ax.set_xlim(xlim)
ax.set_ylim()
plt.show()