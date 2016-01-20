import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/Users/poulin/Documents/Labo/ProgrammeDarkAges/class_public-2.4.2/output_decay/EMDecay_Tau1e25_f1_1GeV_thermodynamics.dat', '/Users/poulin/Documents/Labo/ProgrammeDarkAges/class_public-2.4.2/output_decay/EMDecay_Tau1e25_f1_1TeV_thermodynamics.dat', '/Users/poulin/Documents/Labo/ProgrammeDarkAges/class_public-2.4.2/output_decay/EMDecay_Tau1e25_f1_1GeV_thermodynamics.dat', '/Users/poulin/Documents/Labo/ProgrammeDarkAges/class_public-2.4.2/output_decay/EMDecay_Tau1e26_f1_thermodynamics.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
roots = ['EMDecay_Tau1e25_f1_1GeV_thermodynamics', 'EMDecay_Tau1e25_f1_1TeV_thermodynamics', 'EMDecay_Tau1e25_f1_1GeV_thermodynamics', 'EMDecay_Tau1e26_f1_thermodynamics']

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

index, curve = 2, data[2]
y_axis = [u'Tb[K]']
tex_names = ['Tb [K]']
x_axis = 'z'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 6]))

index, curve = 3, data[3]
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