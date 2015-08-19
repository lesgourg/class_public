import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/Users/poulin/Documents/Doc_Labo/ProgrammeDarkAges/class_public-2.4.2/output_test/Reio_DM_steplike_cl.dat', '/Users/poulin/Documents/Doc_Labo/ProgrammeDarkAges/class_public-2.4.2/output_test/Reio_None_cl.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
roots = ['Reio_DM_steplike_cl', 'Reio_None_cl']

fig, ax = plt.subplots()
y_axis = [u'TT']
tex_names = ['TT']
x_axis = 'l'
ax.set_xlabel('$\ell$', fontsize=16)
plt.show()