import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/Users/poulin/Documents/Labo/ProgrammeDarkAges/class_public-2.4.2/output_5fev/F1_G100_pk.dat', '/Users/poulin/Documents/Labo/ProgrammeDarkAges/class_public-2.4.2/output_5fev/F1_G100_check_pk.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
roots = ['F1_G100_pk', 'F1_G100_check_pk']

fig, ax = plt.subplots()
y_axis = [u'P(Mpc/h)^3']
tex_names = ['P (Mpc/h)^3']
x_axis = 'k (h/Mpc)'
ax.set_xlabel('k (h/Mpc)', fontsize=16)
plt.show()