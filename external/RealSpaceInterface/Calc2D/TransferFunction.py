import os.path
import pickle
import uuid
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline, RectBivariateSpline
import sys
import logging

from classy import Class

import Calc2D.Database as Database
import config

TRANSFER_QUANTITIES = ["d_g", "d_ur", "d_cdm", "d_b", "d_g/4 + psi"]

def ComputeTransferData(settings, redshift):
    database_key = settings.copy()
    database_key.update({'redshift': tuple(redshift)})

    database = Database.Database(config.DATABASE_DIR)
    if database_key in database:
        return database[database_key], redshift
    else:
        cosmo = Class()
        cosmo.set(settings)
        cosmo.compute()

        outputData = [cosmo.get_transfer(z) for z in redshift]
        # Calculate d_g/4+psi
        for transfer_function_dict in outputData:
            transfer_function_dict["d_g/4 + psi"] = transfer_function_dict["d_g"]/4 + transfer_function_dict["psi"]
        # Now filter the relevant fields
        fields = TRANSFER_QUANTITIES + ["k (h/Mpc)"]
        outputData = [{field: outputData[i][field] for field in fields} for i in range(len(redshift))]

        database[database_key] = outputData
        return outputData, redshift


def ComputeTransferFunctionList(cosmologicalParameters, redshift, kperdecade=200, P_k_max=100):
    class_settings = cosmologicalParameters.copy()
    class_settings.update({
        "output": "mTk",
        "gauge": "newtonian",
        "matter_source_in_current_gauge": "yes",
        "evolver": "1",
        "P_k_max_h/Mpc": P_k_max,
        "k_per_decade_for_pk": kperdecade,
        "z_max_pk": str(max(redshift)),
    })

    data_dict, redshift = ComputeTransferData(class_settings, redshift)
    transfer_functions = {field: [] for field in TRANSFER_QUANTITIES}


    for i in range(len(redshift)):
        k_data = data_dict[0]["k (h/Mpc)"] * cosmologicalParameters["h"]  #in order to get k [1/Mpc]
        k_data_zero = np.concatenate(([0.0], k_data))
        for field in TRANSFER_QUANTITIES:
            data = data_dict[i][field] / data_dict[i][field][0]
            data_zero = np.concatenate(([1.0], data))
            interpolated_func = InterpolatedUnivariateSpline(k_data_zero, data_zero)
            transfer_functions[field].append(interpolated_func)

    return transfer_functions
