import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/Users/poulin/Documents/Labo/ProgrammeDarkAges/class_public-2.4.2/test/standard_cl.dat', '/Users/poulin/Documents/Labo/ProgrammeDarkAges/class_public-2.4.2/test/test_halo_master_cl.dat', '/Users/poulin/Documents/Labo/ProgrammeDarkAges/class_public-2.4.2/test/test_halo_master_beyond_cl.dat', '/Users/poulin/Documents/Labo/ProgrammeDarkAges/class_public-2.4.2/test/test_halo_perso_cl.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
roots = ['standard_cl', 'test_halo_master_cl', 'test_halo_master_beyond_cl', 'test_halo_perso_cl']

fig, ax = plt.subplots()
y_axis = [u'TT']
tex_names = ['TT']
x_axis = 'l'
y_axis = [u'TT']
tex_names = ['TT']
x_axis = 'l'
y_axis = [u'TT']
tex_names = ['TT']
x_axis = 'l'
ax.set_xlabel('$\ell$', fontsize=16)
plt.show()