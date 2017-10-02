import matplotlib.pyplot as plt
import numpy as np
import sys
from scipy.interpolate import InterpolatedUnivariateSpline

fig, ax = plt.subplots()
ax.set_xscale('log')
ax.set_yscale('log')


if(len(sys.argv) != 3 and len(sys.argv) != 5):
  print "invalid number of arguments. Use:"
  print "lcdm_pk.dat ncdm_pk.dat"
  print "or"
  print "lcdm_pk.dat alpha beta gamma"
  sys.exit(0)

kmin = 0.5
kmax = 20

def T(a,b,c,x):
  return (1 + (a*x)**b)**c

Plcdm = np.loadtxt(sys.argv[1])
kPlcdm = InterpolatedUnivariateSpline(Plcdm[:,0],Plcdm[:,0]*Plcdm[:,1],k=5)

if(len(sys.argv) == 3):
  P = np.loadtxt(sys.argv[2])
  k = P[:,0]
  kP = InterpolatedUnivariateSpline(P[:,0],P[:,0]*P[:,1],k=5)
else:
  k = Plcdm[:,0]
  P = T(float(sys.argv[2]),float(sys.argv[3]),float(sys.argv[4]),k)**2*Plcdm[:,1]
  kP = InterpolatedUnivariateSpline(k,k*P,k=5)

P1Dlcdm = []
P1D = []
for i in k[:-1]:
  P1Dlcdm.append(1/(2*np.pi)*kPlcdm.integral(i,'inf'))
  P1D.append(1/(2*np.pi)*kP.integral(i,'inf'))
P1Dlcdm = np.array(P1Dlcdm)
P1D = np.array(P1D)

ax.plot(k[:-1],np.abs(P1Dlcdm))
ax.plot(k[:-1],np.abs(P1D))

r = P1D/P1Dlcdm
ax.plot(k[:-1],r)
r = InterpolatedUnivariateSpline(k[:-1],r,k=5)

A = r.integral(kmin, kmax)
Alcdm = kmax - kmin

dA = (Alcdm - A)/Alcdm

print dA

#plt.show()
