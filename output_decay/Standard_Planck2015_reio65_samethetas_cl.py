import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/Users/poulin/Documents/Labo/ProgrammeDarkAges/class_public-2.4.2/output_decay/Standard_Planck2015_reio65_samethetas_cl.dat', '/Users/poulin/Documents/Labo/ProgrammeDarkAges/class_public-2.4.2/output_decay/EMDecay_Tau1e15_f1emoins6_D1_samethetas_cl.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
roots = ['Standard_Planck2015_reio65_samethetas_cl', 'EMDecay_Tau1e15_f1emoins6_D1_samethetas_cl']

fig, ax = plt.subplots()
y_axis = [u'EE']
tex_names = ['EE']
x_axis = 'l'
ax.set_xlabel('$\ell$', fontsize=16)
plt.show()