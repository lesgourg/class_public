import numpy as np
import numpy.fft as fft

def realFourier(step, Value):
    FValue = np.fft.fftshift(
        np.fft.rfft2(Value), axes=(0))  #shifting only the x axes

    kx = np.fft.fftshift(np.fft.fftfreq(Value.shape[0], d=step)) * 2 * np.pi
    ky = np.fft.rfftfreq(Value.shape[0], d=step) * 2 * np.pi

    return kx, ky, FValue

def realInverseFourier(FValue):
    return np.fft.irfft2(np.fft.ifftshift(
        FValue, axes=(0)))  #shifting only on the x axes


def realInverseAllFourier(allFValue):
    return np.fft.irfftn(
        np.fft.ifftshift(allFValue, axes=(1)),
        axes=(1, 2))  #shifting only on the x axes
