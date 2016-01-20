import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/Users/poulin/Documents/Labo/ProgrammeDarkAges/class_public-2.4.2/output_decay/EMDecay_Tau1e15_f1emoins7_D1_samethetas_thermodynamics.dat', '/Users/poulin/Documents/Labo/ProgrammeDarkAges/class_public-2.4.2/output_decay/EMDecay_Tau1e15_f1emoins7_D1_withT_thermodynamics.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
roots = ['EMDecay_Tau1e15_f1emoins7_D1_samethetas_thermodynamics', 'EMDecay_Tau1e15_f1emoins7_D1_withT_thermodynamics']

fig, ax = plt.subplots()

index, curve = 0, data[0]
y_axis = [u'Tb[K]']
tex_names = ['Tb [K]']
x_axis = 'z'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 6]))

index, curve = 1, data[1]
y_axis = [u'Tb[K]']
tex_names = ['Tb [K]']
x_axis = 'z'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 6]))

ax.legend([root+': '+elem for (root, elem) in
    itertools.product(roots, y_axis)], loc='best')

ax.set_xlabel('z', fontsize=16)
plt.show()