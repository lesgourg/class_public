import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/Users/poulin/Documents/Doc_Labo/ProgrammeDarkAges/class_public-2.4.2/output_decay/test_EMDecay_G10_Omega20_cl.dat', '/Users/poulin/Documents/Doc_Labo/ProgrammeDarkAges/class_public-2.4.2/output_decay/EMDecay_G10_D0_cl.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
roots = ['test_EMDecay_G10_Omega20_cl', 'EMDecay_G10_D0_cl']

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

ax.legend([root+': '+elem for (root, elem) in
    itertools.product(roots, y_axis)], loc='best')

ax.set_xlabel('$\ell$', fontsize=16)
plt.show()