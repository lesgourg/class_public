import os
import logging

import cv2
import numpy as np

from classy import Class

from Calc2D.TransferFunction import ComputeTransferFunctionList
from Calc2D.DataGeneration import GenerateGaussianData, GenerateSIData
from Calc2D.DataPropagation import PropagateDatawithList
from Calc2D.rFourier import *
from Calc2D.Database import Database
from collections import namedtuple

import config

ClSpectrum = namedtuple("Cls", ["l", "tCl"])
PkSpectrum = namedtuple("Pkh", ["kh", "Pkh"])

def normalize(real):
    """
    Given the `real` data, i.e. either a 2d array or a flattened 1d array
    of the relative density perturbations, normalize its values as follows:

    The client expects the values to be in the interval of [-1, 1].
    Take a symmetric interval around a `real` value of 0 and linearly
    map it to the required interval [-1, 1].
    """
    minimum, maximum = real.min(), real.max()
    bound = max(abs(maximum), abs(minimum))
    result = real / bound
    return result

class Calculation(object):
    def __init__(self,
                 kbins,
                 resolution = 200,
                 gauge="newtonian",
                 kperdecade=200,
                 P_k_max=100,
                 evolver=1):
        # also sets `endshape` through setter
        self.resolution = resolution
        self.gauge = gauge
        self.evolver = evolver
        self.P_k_max = P_k_max
        self.kperdecade = kperdecade

        self.redshift = None # to be set later
        self.size = None # to be set later

        self.krange = np.logspace(-4, 1, kbins)

    @property
    def resolution(self):
        return self._resolution

    @resolution.setter
    def resolution(self, resolution):
        self._resolution = resolution
        self.endshape = (resolution, resolution)

    def getData(self, redshiftindex):
        FValuenew = PropagateDatawithList(
            k=self.k,
            FValue=self.FValue,
            zredindex=redshiftindex,
            transferFunctionlist=self.TransferFunctionList)

        Valuenew = dict()
        FValue_abs = np.abs(self.FValue)
        _min, _max = FValue_abs.min(), FValue_abs.max()
        dimensions = (self.endshape[0] / 2, self.endshape[1])
        for quantity, FT in FValuenew.items():
            FT_abs = np.abs(FT)
            FT_normalized = cv2.resize(FT_abs, dimensions).ravel()
            FT_normalized = (FT_normalized - _min) / (_max - _min)
            real = realInverseFourier(FT.reshape(self.FValue.shape))
            # real = cv2.resize(real, self.endshape).ravel()
            real = real.ravel()

            minimum, maximum = real.min(), real.max()
            Valuenew[quantity] = normalize(real)

        return Valuenew, FValuenew, (minimum, maximum)


    def getInitialData(self):
        # for odd values of self._resolution, this is necessary
        Value = cv2.resize(realInverseFourier(self.FValue), self.endshape)

        minimum, maximum = Value.min(), Value.max()
        Value = normalize(Value)

        assert Value.size == self.resolution ** 2

        return Value.ravel(), cv2.resize(
            (np.abs(self.FValue) - np.abs(self.FValue).min()) /
            (np.abs(self.FValue).max() - np.abs(self.FValue).min()),
            (self.endshape[0] / 2, self.endshape[1])).ravel(), (minimum,
                                                                maximum)

    def getTransferData(self, redshiftindex):
        return {field: transfer_function[redshiftindex](self.krange) for field, transfer_function in self.TransferFunctionList.items()}, self.krange

    def setCosmologialParameters(self, cosmologicalParameters):
        self.cosmologicalParameters = cosmologicalParameters

        # Calculate transfer functions
        self.TransferFunctionList = ComputeTransferFunctionList(self.cosmologicalParameters, self.redshift)
        # Calculate Cl's
        self.tCl, self.mPk = self.calculate_spectra(self.cosmologicalParameters)

    def calculate_spectra(self, cosmo_params, force_recalc=False):
        settings = cosmo_params.copy()
        settings.update({
            "output": "tCl,mPk",
            "evolver": "1",
            "gauge": "newtonian",
            "P_k_max_1/Mpc": 10,
            })

        database = Database(config.DATABASE_DIR, "spectra.dat")

        if settings in database and not force_recalc:
            data = database[settings]
            ell = data["ell"]
            tt = data["tt"]
            kh = data["kh"]
            Pkh = data["Pkh"]
            self.z_rec = data["z_rec"]
        else:
            cosmo = Class()
            cosmo.set(settings)
            cosmo.compute()
            # Cl's
            data = cosmo.raw_cl()
            ell = data["ell"]
            tt = data["tt"]
            # Matter spectrum
            k = np.logspace(-3, 1, config.MATTER_SPECTRUM_CLIENT_SAMPLES_PER_DECADE * 4)
            Pk = np.vectorize(cosmo.pk)(k, 0)
            kh = k * cosmo.h()
            Pkh = Pk / cosmo.h()**3
            # Get redshift of decoupling
            z_rec = cosmo.get_current_derived_parameters(['z_rec'])['z_rec']
            self.z_rec = z_rec
            # Store to database
            database[settings] = {
            "ell": data["ell"],
            "tt": data["tt"],

            "kh": k,
            "Pkh": Pk,

            "z_rec": z_rec,
            }

        return ClSpectrum(ell[2:], tt[2:]), PkSpectrum(kh, Pkh)

    @property
    def z_dec(self):
        if self.z_rec is None:
            raise ValueError("z_rec hasn't been computed yet")
        return self.z_rec

    def setInitialConditions(self,
                             A=1,
                             sigma=2,
                             initialDataType="SI",
                             SIlimit=None,
                             SI_ns=0.96):
        logging.info("Generating Initial Condition")

        if initialDataType == "Gaussian":
            self.ValueE, self.FValue, self.k, self.kxE, self.kyE = GenerateGaussianData(
                sigma, self.size, self.resolution)
        elif initialDataType == "SI":
            self.ValueE, self.FValue, self.k, self.kxE, self.kyE = GenerateSIData(
                A,
                self.size,
                self.resolution,
                limit=SIlimit,
                ns=SI_ns)
        else:
            logging.warn("initialDataType " + str(initialDataType) + " not found")
