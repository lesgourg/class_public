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
  return (1. + (a*x)**b)**c

Plcdm_file = np.loadtxt(sys.argv[1])
k=Plcdm_file[:,0]
Plcdm=Plcdm_file[:,1]
kPlcdm = InterpolatedUnivariateSpline(Plcdm_file[:,0],Plcdm_file[:,0]*Plcdm_file[:,1],k=5)

if(len(sys.argv) == 3):
  P_file = np.loadtxt(sys.argv[2])
  interp=InterpolatedUnivariateSpline(P_file[:,0],P_file[:,1],k=5)
  P=interp(k)
  kP = InterpolatedUnivariateSpline(k,k*P,k=5)
else:
  P = T(float(sys.argv[2]),float(sys.argv[3]),float(sys.argv[4]),k)**2*Plcdm_file[:,1]
  kP = InterpolatedUnivariateSpline(k,k*P,k=5)

Pratio=P/Plcdm
Phalf=Pratio[np.abs(Pratio-1./2.).argmin()]
ihalf=np.argwhere(Pratio==Phalf)
k_half=np.asscalar(k[ihalf])
print 'k_1/2',k_half

P1Dlcdm = []
P1D = []
for i in k[:-1]:
  P1Dlcdm.append(1./(2.*np.pi)*kPlcdm.integral(i,'inf'))
  P1D.append(1./(2.*np.pi)*kP.integral(i,'inf'))
P1Dlcdm = np.array(P1Dlcdm)
P1D = np.array(P1D)

ax.plot(k[:-1],np.abs(P1Dlcdm))
ax.plot(k[:-1],np.abs(P1D))

r = P1D/P1Dlcdm
#ax.plot(k[:-1],r)
r = InterpolatedUnivariateSpline(k[:-1],r,k=5)
k_int=np.linspace(0.5,20.,1000)
r2 = r(k_int)
r_int = InterpolatedUnivariateSpline(k_int,r2,k=5)

A = r_int.integral(kmin, kmax)
Alcdm = kmax - kmin

dA = (Alcdm - A)/Alcdm

print 'dA',dA

ax.set_xscale('log')
ax.set_xlim([0.1,20])
ax.set_yscale('log')
ax.set_ylim([1.e-4,1.e4])
plt.show()
plt.savefig('plot_lyman.pdf')
