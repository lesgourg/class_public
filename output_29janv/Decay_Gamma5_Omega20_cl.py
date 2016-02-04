import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/Users/poulin/Documents/Labo/ProgrammeDarkAges/class_public-2.4.2/output_29janv/Decay_Gamma5_Omega20_cl.dat', '/Users/poulin/Documents/Labo/ProgrammeDarkAges/class_public-2.4.2/output_29janv/Decay_Gamma1_Omega100_cl.dat', '/Users/poulin/Documents/Labo/ProgrammeDarkAges/class_public-2.4.2/output_29janv/Decay_Gamma25_Omega20_cl.dat', '/Users/poulin/Documents/Labo/ProgrammeDarkAges/class_public-2.4.2/output_29janv/Decay_Gamma5_Omega100_cl.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
roots = ['Decay_Gamma5_Omega20_cl', 'Decay_Gamma1_Omega100_cl', 'Decay_Gamma25_Omega20_cl', 'Decay_Gamma5_Omega100_cl']

fig, ax = plt.subplots()

index, curve = 0, data[0]
y_axis = [u'TT']
tex_names = ['TT']
x_axis = 'l'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 1]))

index, curve = 1, data[1]
y_axis = [u'TT']
tex_names = ['TT']
x_axis = 'l'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 1]))

index, curve = 2, data[2]
y_axis = [u'TT']
tex_names = ['TT']
x_axis = 'l'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 1]))

index, curve = 3, data[3]
y_axis = [u'TT']
tex_names = ['TT']
x_axis = 'l'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 1]))

ax.legend([root+': '+elem for (root, elem) in
    itertools.product(roots, y_axis)], loc='best')

ax.set_xlabel('$\ell$', fontsize=16)
plt.show()