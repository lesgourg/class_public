import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/Users/poulin/Documents/Doc_Labo/ProgrammeDarkAges/class_public-2.4.2/output_new/electrons_1GeV_zh10_B4_taureio79_thermodynamics.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
roots = ['electrons_1GeV_zh10_B4_taureio79_thermodynamics']

fig, ax = plt.subplots()

index, curve = 0, data[0]
y_axis = [u'exp(-kappa)']
tex_names = ['exp(-kappa)']
x_axis = 'z'
ylim = [0.0001, 1.0]
xlim = [1.0, 600.0]
ax.loglog(curve[:, 0], abs(curve[:, 4]))

ax.legend([root+': '+elem for (root, elem) in
    itertools.product(roots, y_axis)], loc='best')

ax.set_xlabel('z', fontsize=16)
ax.set_xlim(xlim)
ax.set_ylim()
ax.set_ylim(ylim)
plt.show()