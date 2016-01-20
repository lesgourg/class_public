import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/Users/poulin/Documents/Labo/ProgrammeDarkAges/class_public-2.4.2/output_decay/Standard_Planck2015_reio65_thermodynamics.dat', '/Users/poulin/Documents/Labo/ProgrammeDarkAges/class_public-2.4.2/output_decay/EMDecay_Tau1e15_f1emoins7_D1_beyond_thermodynamics.dat', '/Users/poulin/Documents/Labo/ProgrammeDarkAges/class_public-2.4.2/output_decay/EMDecay_Tau1e15_f1emoins7_D1_onthespot_thermodynamics.dat', '/Users/poulin/Documents/Labo/ProgrammeDarkAges/class_public-2.4.2/output_decay/EMDecay_Tau1e25_f1_D1_onthespot_thermodynamics.dat', '/Users/poulin/Documents/Labo/ProgrammeDarkAges/class_public-2.4.2/output_decay/EMDecay_Tau1e25_f1_D1_beyond_thermodynamics.dat', '/Users/poulin/Documents/Labo/ProgrammeDarkAges/class_public-2.4.2/output_decay/EMDecay_Tau1e25_f1emoins7_D1_beyond1TeV_thermodynamics.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
roots = ['Standard_Planck2015_reio65_thermodynamics', 'EMDecay_Tau1e15_f1emoins7_D1_beyond_thermodynamics', 'EMDecay_Tau1e15_f1emoins7_D1_onthespot_thermodynamics', 'EMDecay_Tau1e25_f1_D1_onthespot_thermodynamics', 'EMDecay_Tau1e25_f1_D1_beyond_thermodynamics', 'EMDecay_Tau1e25_f1emoins7_D1_beyond1TeV_thermodynamics']

fig, ax = plt.subplots()

index, curve = 0, data[0]
y_axis = [u'x_e']
tex_names = ['x_e']
x_axis = 'z'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 2]))

index, curve = 1, data[1]
y_axis = [u'x_e']
tex_names = ['x_e']
x_axis = 'z'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 2]))

index, curve = 2, data[2]
y_axis = [u'x_e']
tex_names = ['x_e']
x_axis = 'z'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 2]))

index, curve = 3, data[3]
y_axis = [u'x_e']
tex_names = ['x_e']
x_axis = 'z'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 2]))

index, curve = 4, data[4]
y_axis = [u'x_e']
tex_names = ['x_e']
x_axis = 'z'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 2]))

index, curve = 5, data[5]
y_axis = [u'x_e']
tex_names = ['x_e']
x_axis = 'z'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 2]))

ax.legend([root+': '+elem for (root, elem) in
    itertools.product(roots, y_axis)], loc='best')

ax.set_xlabel('z', fontsize=16)
plt.show()