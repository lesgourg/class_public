import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/Users/poulin/Documents/Labo/ProgrammeDarkAges/class_public-2.4.2/output_4fev/check_decay_noDCDM_pk.dat', '/Users/poulin/Documents/Labo/ProgrammeDarkAges/class_public-2.4.2/output_4fev/check_decay_pk.dat', '/Users/poulin/Documents/Labo/ProgrammeDarkAges/class_public-2.4.2/output_29janv/Standard_Planck2015_pk.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
roots = ['check_decay_noDCDM_pk', 'check_decay_pk', 'Standard_Planck2015_pk']

fig, ax = plt.subplots()

index, curve = 0, data[0]
y_axis = [u'P(Mpc/h)^3']
tex_names = ['P (Mpc/h)^3']
x_axis = 'k (h/Mpc)'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 1]))

index, curve = 1, data[1]
y_axis = [u'P(Mpc/h)^3']
tex_names = ['P (Mpc/h)^3']
x_axis = 'k (h/Mpc)'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 1]))

index, curve = 2, data[2]
y_axis = [u'P(Mpc/h)^3']
tex_names = ['P (Mpc/h)^3']
x_axis = 'k (h/Mpc)'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 1]))

ax.legend([root+': '+elem for (root, elem) in
    itertools.product(roots, y_axis)], loc='best')

ax.set_xlabel('k (h/Mpc)', fontsize=16)
plt.show()