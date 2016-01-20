import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/Users/poulin/Documents/Labo/ProgrammeDarkAges/class_public-2.4.2/output_decay/EMDecay_Tau1e14_f1emoins8_D1_onthespot_thermodynamics.dat', '/Users/poulin/Documents/Labo/ProgrammeDarkAges/class_public-2.4.2/output_decay/EMDecay_Tau1e15_f1emoins7_D1_onthespot_thermodynamics.dat', '/Users/poulin/Documents/Labo/ProgrammeDarkAges/class_public-2.4.2/output_decay/EMDecay_Tau1e13_f1emoins8_D1_onthespot_thermodynamics.dat', '/Users/poulin/Documents/Labo/ProgrammeDarkAges/class_public-2.4.2/output_decay/EMDecay_Tau1e15_f1emoins6_D1_onthespot_thermodynamics.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
roots = ['EMDecay_Tau1e14_f1emoins8_D1_onthespot_thermodynamics', 'EMDecay_Tau1e15_f1emoins7_D1_onthespot_thermodynamics', 'EMDecay_Tau1e13_f1emoins8_D1_onthespot_thermodynamics', 'EMDecay_Tau1e15_f1emoins6_D1_onthespot_thermodynamics']

fig, ax = plt.subplots()

index, curve = 0, data[0]
y_axis = []
tex_names = []
x_axis = 'z'
ylim = []
xlim = []

index, curve = 1, data[1]
y_axis = []
tex_names = []
x_axis = 'z'
ylim = []
xlim = []

index, curve = 2, data[2]
y_axis = []
tex_names = []
x_axis = 'z'
ylim = []
xlim = []

index, curve = 3, data[3]
y_axis = []
tex_names = []
x_axis = 'z'
ylim = []
xlim = []

ax.legend([root+': '+elem for (root, elem) in
    itertools.product(roots, y_axis)], loc='best')

ax.set_xlabel('z', fontsize=16)
plt.show()