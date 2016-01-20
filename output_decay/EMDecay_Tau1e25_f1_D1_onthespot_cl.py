import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/Users/poulin/Documents/Labo/ProgrammeDarkAges/class_public-2.4.2/output_decay/EMDecay_Tau1e25_f1_D1_onthespot_cl.dat', '/Users/poulin/Documents/Labo/ProgrammeDarkAges/class_public-2.4.2/output_decay/EMDecay_Tau1e25_f1_D1_beyond_cl.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
roots = ['EMDecay_Tau1e25_f1_D1_onthespot_cl', 'EMDecay_Tau1e25_f1_D1_beyond_cl']

fig, ax = plt.subplots()
y_axis = [u'TT']
tex_names = ['TT']
x_axis = 'l'
ax.set_xlabel('$\ell$', fontsize=16)
plt.show()