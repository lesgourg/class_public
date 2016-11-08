#!/usr/bin/python
from __future__ import print_function
import sys
from math import exp

# README:
#
# This is an example python script for the external_Pk mode of Class.
# It generates the primordial spectrum of LambdaCDM.
# It can be edited and used directly, though keeping a copy of it is recommended.
#
# Two (maybe three) things need to be edited:
#
# 1. The name of the parameters needed for the calculation of Pk.
#    "sys.argv[1]" corresponds to "custom1" in Class, an so on

try :
    k_0           = float(sys.argv[1])
    A             = float(sys.argv[2])
    n_s           = float(sys.argv[3])

# Error control, no need to touch
except IndexError :
    raise IndexError("It seems you are calling this script with too few arguments.")
except ValueError :
    raise ValueError("It seems some of the arguments are not correctly formatted. "+
                     "Remember that they must be floating point numbers.")

# 2. The function giving P(k), including the necessary import statements.
#    Inside this function, you can use the parameters named in the previous step.

def P(k) :
    return A * (k/k_0)**(n_s-1.)

# 3. Limits for k and precision:
#    Check that the boundaries are correct for your case.
#    It is safer to set k_per_decade primordial slightly bigger than that of Class.

k_min  = 1.e-6
k_max  = 10.
k_per_decade_primordial = 200.

#
# And nothing should need to be edited from here on.
#

# Filling the array of k's
ks = [float(k_min)]
while ks[-1] <= float(k_max) :
    ks.append(ks[-1]*10.**(1./float(k_per_decade_primordial)))

# Filling the array of Pk's
for k in ks :
    P_k = P(k)
    print("%.18g %.18g" % (k, P_k))

