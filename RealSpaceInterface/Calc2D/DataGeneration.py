import logging

import numpy as np
import cv2

from Calc2D.rFourier import realFourier, realInverseFourier

def GenerateGaussianData(sigma, size, points, A=1):
    xr = np.linspace(-size / 2.0, size / 2.0, points)
    yr = np.linspace(-size / 2.0, size / 2.0, points)
    step = xr[1] - xr[0]
    x, y = np.meshgrid(
        xr, yr, indexing='ij', sparse=True)  # indexing is important
    del xr, yr

    #use the more easy formula
    Value = A * np.exp(-(x**2 + y**2) / (2 * sigma**2))

    kx, ky, FValue = realFourier(step, Value)
    kxr, kyr = np.meshgrid(kx, ky, indexing='ij', sparse=True)

    k = np.sqrt(kxr**2 + kyr**2)
    del kxr, kyr

    kx = (min(kx), max(kx))  #just return the extremal values to save memory
    ky = (min(ky), max(ky))

    ValueE = (Value.min(), Value.max())

    return ValueE, FValue, k, kx, ky

def GenerateSIData(A, size, points, limit=None, ns=0.96):
    xr = np.linspace(-size / 2.0, size / 2.0, points)
    yr = np.linspace(-size / 2.0, size / 2.0, points)
    step = xr[1] - xr[0]

    x, y = np.meshgrid(
        xr, yr, indexing='ij', sparse=True)  # indexing is important
    del xr, yr
    Value = 0 * x + 0 * y

    kx, ky, FValue = realFourier(step, Value)  #FValue==0

    kxr, kyr = np.meshgrid(kx, ky, indexing='ij', sparse=True)

    k = np.sqrt(kxr**2 + kyr**2)
    del kxr, kyr

    if limit == None:

        ktilde = k.flatten()
        ktilde[np.argmin(k)] = 10**9  #just let the background be arbitrary low
        ktilde = ktilde.reshape(k.shape)

        FValue = np.random.normal(
            loc=0,
            scale=np.sqrt(A / ktilde**(
                2 - (ns - 1) * 2. / 3.)) / np.sqrt(2)) + np.random.normal(
                    loc=0,
                    scale=np.sqrt(A / ktilde**
                                  (2 - (ns - 1) * 2. / 3.)) / np.sqrt(2)) * 1j

    elif type(limit) == list or type(limit) == tuple:

        iunder, junder = np.where(k < limit[1])

        for t in range(len(iunder)):

            if k[iunder[t]][junder[t]] > limit[0] and k[iunder[t]][junder[t]] > 0:

                FValue[iunder[t]][junder[t]] = np.random.normal(
                    loc=0,
                    scale=np.sqrt(A / k[iunder[t]][junder[t]]**
                                  (2 - (ns - 1) * 2. / 3.)) /
                    np.sqrt(2)) + np.random.normal(
                        loc=0,
                        scale=np.sqrt(A / k[iunder[t]][junder[t]]**
                                      (2 -
                                       (ns - 1) * 2. / 3.)) / np.sqrt(2)) * 1j

    else:
        raise ValueError("limit must be None or tuple or list")

    Value = realInverseFourier(FValue)

    kx = (min(kx), max(kx))
    ky = (min(ky), max(ky))

    ValueE = (Value.min(), Value.max())

    return ValueE, FValue, k, kx, ky
