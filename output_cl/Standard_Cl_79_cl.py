import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/Users/poulin/Documents/Doc_Labo/ProgrammeDarkAges/class_public-2.4.2/output_cl/Standard_Cl_79_cl.dat', '/Users/poulin/Documents/Doc_Labo/ProgrammeDarkAges/class_public-2.4.2/output_cl/Muons_1GeV_NoHalo_0079_cl.dat', '/Users/poulin/Documents/Doc_Labo/ProgrammeDarkAges/class_public-2.4.2/output_cl/Muons_1GeV_zh30_B8_0079_cl.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
roots = ['Standard_Cl_79_cl', 'Muons_1GeV_NoHalo_0079_cl', 'Muons_1GeV_zh30_B8_0079_cl']

fig, ax = plt.subplots()
y_axis = [u'TT']
tex_names = ['TT']
x_axis = 'l'
y_axis = [u'TT']
tex_names = ['TT']
x_axis = 'l'
ax.set_xlabel('$\ell$', fontsize=16)
plt.show()