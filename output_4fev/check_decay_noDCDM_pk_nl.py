import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/Users/poulin/Documents/Labo/ProgrammeDarkAges/class_public-2.4.2/output_4fev/check_decay_noDCDM_pk_nl.dat', '/Users/poulin/Documents/Labo/ProgrammeDarkAges/class_public-2.4.2/output_4fev/check_decay_pk_nl.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
roots = ['check_decay_noDCDM_pk_nl', 'check_decay_pk_nl']

fig, ax = plt.subplots()
y_axis = [u'P(Mpc/h)^3']
tex_names = ['P (Mpc/h)^3']
x_axis = 'k (h/Mpc)'
ax.set_xlabel('k (h/Mpc)', fontsize=16)
plt.show()