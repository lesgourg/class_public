import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/Users/poulin/Documents/Labo/ProgrammeDarkAges/class_public-2.4.2/output_29janv/Standard_Planck2015_cl.dat', '/Users/poulin/Documents/Labo/ProgrammeDarkAges/class_public-2.4.2/output_5fev/EMDecay_ee1GeV_tau1e25_onthespot_hyrecv2_cl.dat', '/Users/poulin/Documents/Labo/ProgrammeDarkAges/class_public-2.4.2/output_5fev/EMDecay_f1_tau1e25_beyond_hyrecv2_cl.dat', '/Users/poulin/Documents/Labo/ProgrammeDarkAges/class_public-2.4.2/output_5fev/EMDecay_ee1GeV_tau1e15_onthespot_hyrecv2_cl.dat', '/Users/poulin/Documents/Labo/ProgrammeDarkAges/class_public-2.4.2/output_5fev/EMDecay_ee1GeV_tau1e15_beyond_hyrecv2_cl.dat', '/Users/poulin/Documents/Labo/ProgrammeDarkAges/class_public-2.4.2/output_5fev/EMDecay_fmoins6_tau1e18_onthespot_hyrecv2_cl.dat', '/Users/poulin/Documents/Labo/ProgrammeDarkAges/class_public-2.4.2/output_5fev/EMDecay_fmoins6_tau1e18_beyond_hyrecv2_cl.dat', '/Users/poulin/Documents/Labo/ProgrammeDarkAges/class_public-2.4.2/output_5fev/EMDecay_fmoins8_tau1e13_onthespot_hyrecv2_cl.dat', '/Users/poulin/Documents/Labo/ProgrammeDarkAges/class_public-2.4.2/output_5fev/EMDecay_fmoins8_tau1e13_beyond_hyrecv2_cl.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
roots = ['Standard_Planck2015_cl', 'EMDecay_ee1GeV_tau1e25_onthespot_hyrecv2_cl', 'EMDecay_f1_tau1e25_beyond_hyrecv2_cl', 'EMDecay_ee1GeV_tau1e15_onthespot_hyrecv2_cl', 'EMDecay_ee1GeV_tau1e15_beyond_hyrecv2_cl', 'EMDecay_fmoins6_tau1e18_onthespot_hyrecv2_cl', 'EMDecay_fmoins6_tau1e18_beyond_hyrecv2_cl', 'EMDecay_fmoins8_tau1e13_onthespot_hyrecv2_cl', 'EMDecay_fmoins8_tau1e13_beyond_hyrecv2_cl']

fig, ax = plt.subplots()

index, curve = 0, data[0]
y_axis = [u'EE']
tex_names = ['EE']
x_axis = 'l'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 2]))

index, curve = 1, data[1]
y_axis = [u'EE']
tex_names = ['EE']
x_axis = 'l'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 2]))

index, curve = 2, data[2]
y_axis = [u'EE']
tex_names = ['EE']
x_axis = 'l'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 2]))

index, curve = 3, data[3]
y_axis = [u'EE']
tex_names = ['EE']
x_axis = 'l'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 2]))

index, curve = 4, data[4]
y_axis = [u'EE']
tex_names = ['EE']
x_axis = 'l'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 2]))

index, curve = 5, data[5]
y_axis = [u'EE']
tex_names = ['EE']
x_axis = 'l'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 2]))

index, curve = 6, data[6]
y_axis = [u'EE']
tex_names = ['EE']
x_axis = 'l'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 2]))

index, curve = 7, data[7]
y_axis = [u'EE']
tex_names = ['EE']
x_axis = 'l'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 2]))

index, curve = 8, data[8]
y_axis = [u'EE']
tex_names = ['EE']
x_axis = 'l'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 2]))

ax.legend([root+': '+elem for (root, elem) in
    itertools.product(roots, y_axis)], loc='best')

ax.set_xlabel('$\ell$', fontsize=16)
plt.show()