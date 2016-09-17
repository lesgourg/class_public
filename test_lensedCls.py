import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/Users/poulin/Documents/Labo/ProgrammeDarkAges/class_public-2.4.2/test_lensedCls.dat', '/Users/poulin/Documents/Labo/ProgrammeDarkAges/class_public-2.4.2/test_lensedCls_T25.dat', '/Users/poulin/Documents/Labo/ProgrammeDarkAges/class_public-2.4.2/test_lensedCls_T3.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
roots = ['test_lensedCls', 'test_lensedCls_T25', 'test_lensedCls_T3']

fig, ax = plt.subplots()
y_axis = [u'TT']
tex_names = ['TT']
x_axis = 'l'
y_axis = [u'TT']
tex_names = ['TT']
x_axis = 'l'
ax.set_xlabel('$\ell$', fontsize=16)
ax.set_xlim(xlim)
ax.set_ylim()
plt.show()