import multiprocessing
from time import time
import sys
import os
import logging
import shutil
import json

import numpy as np
import h5py as h5
from tqdm import tqdm

from classy import Class

from .generate_cosmological_parameters import sample_cosmological_parameters


MiB = 1024**2
GiB = 1024**3
# A (generously) rounded-up estimate of how much space a single
# training file will occupy on disk.
FILE_SIZE_ESTIMATE = 30 * MiB

def generate_data(count, domain, fixed_params, directory, processes=None):
    """
    Generate training/validation/testing set of size `count` by sampling the `domain` of cosmological parameters.

    The files will be placed in the `directory`.

    Optionally, the number of worker processes can be specified using `processes`;
    if not specified, the number of processes will equal the number of CPU cores.

    This function will also save the cosmological parameters to a file "parameters.h5"
    in the same `directory`.

    In addition, the global minima and maxima for each quantity over the
    entire set are written to the file "normalization.json" in `directory`
    as well.

    This function will return a Dict[str, np.ndarray] of the sampled cosmological parameters.
    """
    # Create output directory if it doesn't exist
    os.makedirs(directory, exist_ok=True)
    # Check whether there is enough space available
    expected_size = count * FILE_SIZE_ESTIMATE
    usage = shutil.disk_usage(directory)
    if usage.free < expected_size:
        raise ValueError(
            "You requested a dataset of {} elements which would "
            "require roughly {:.3f} GiB of free space; "
            "however, there are only {:.3f} GiB available!".format(
                count, expected_size / GiB, usage.free / GiB
            ))

    # sample cosmological parameters...
    varying = sample_cosmological_parameters(domain, count)
    # ...and generate training/validation/testing data.
    minima, maxima = generate_source_functions_for(fixed_params, varying, directory, processes=processes)
    # save minima and maxima
    with open(os.path.join(directory, "normalization.json"), "w") as out:
        json.dump({"min": minima, "max": maxima}, out)
    # save cosmological parameters
    with h5.File(os.path.join(directory, "parameters.h5"), "w") as out:
        for key, data in varying.items():
            out.create_dataset(key, data=data)
    # and also return them, just in case
    return varying


def generate_source_functions_for(fixed_params, varying_params, directory, processes=None):
    """
    Generate training/validation/testing data for the given parameters.
    Parameters are specified using the two dicts `fixed_params`, a dict
    containing the _shared_ parameters for all cosmologies, and `varying_params`,
    a dict whose keys are arrays containing the parameters individual to each
    cosmology.
    """
    any_key = next(iter(varying_params))
    count = len(varying_params[any_key])

    # Prepare function arguments
    args = []
    for i in range(count):
        params = fixed_params.copy()
        cosmo_params = {key: varying_params[key][i] for key in varying_params}
        output_file = os.path.join(directory, "sources_{}.h5".format(i))
        arg = (fixed_params, cosmo_params, output_file)
        args.append(arg)

    minima = None
    maxima = None
    with multiprocessing.Pool(processes) as pool:
        for (mins, maxs) in tqdm(pool.imap_unordered(generate_source_function, args), total=count):
            if minima is None:
                minima = mins
                maxima = maxs
            else:
                minima = {key: min(minima[key], mins[key]) for key in mins}
                maxima = {key: max(maxima[key], maxs[key]) for key in maxs}
    return minima, maxima


def generate_source_function(args):
    """
    Generate source function of cosmo_param_sets[i] and stores it into the corresponding file.
    """
    params, cosmo_params, path = args

    logger = logging.getLogger("generate_source_function")

    # create instance of the class "Class"
    cosmo = Class()

    params = params.copy()
    params.update(cosmo_params)
    # TODO don't do this here!
    # Set reasonably high k value
    params["P_k_max_1/Mpc"] = 10.0
    # Make sure that the lowest k value produced by CLASS is _always_ below
    # the lowest k value of our standard k sampling; for this, we set the
    # k_min_tau0 value (the smallest value of k * tau0) to 0.1 of its default
    # value
    params["k_min_tau0"] = 0.01

    cosmo.set(**params)

    # run CLASS
    logger.debug("Computing source functions...", end="")
    cosmo.compute(level=['perturb'])
    logger.debug("done!")

    # get sources from class
    sources_dict, k, tau = cosmo.get_sources()
    source_names = sources_dict.keys()

    # OPTIONAL: Select only a subset of source functions to include
    # in order to reduce file size
    source_names = [
            "t0_isw", "t0_reco_no_isw", "t0_reio_no_isw",
            "t1",
            "t2", "t2_reco", "t2_reio",
            "phi_plus_psi", "delta_m"]

    # create new hdf5 file
    file = h5.File(path, "w")
    file.create_group('sources')
    file.create_group('thermos')
    file.create_group('cosmos')
    file.create_group("sampling")
    file.create_group("background")

    file.create_dataset("sampling/k", data=k)
    file.create_dataset("sampling/tau", data=tau)

    # fill file with source function values
    for name in source_names:
        source = sources_dict[name]
        logger.debug("-- Writing source function", name)
        file.create_dataset("sources/{}".format(name), data=source)

    for name, value in cosmo_params.items():
        file.create_dataset("cosmos/{}".format(name), data=value)

    bg_nn = cosmo.get_backgrounds_for_NN()

    tau_bg = bg_nn["tau"]
    r_s_bg = bg_nn["r_s"]
    rho_b = bg_nn["rho_b"]
    rho_g = bg_nn["rho_g"]

    thermos_NN = cosmo.get_thermos_for_NN()

    thermodynamic_parameters = cosmo.get_current_derived_parameters(['tau_rec',
                                                                        'rs_rec',
                                                                        'ds_rec',
                                                                        'ra_rec',
                                                                        'da_rec',
                                                                        'rd_rec',
                                                                        'z_reio'
                                                                        ])

    file.create_dataset('thermos/tau_reio',data=cosmo.tau_of_z(thermodynamic_parameters['z_reio']))
    file.create_dataset('thermos/tau_rec',data=thermodynamic_parameters['tau_rec'])
    file.create_dataset('thermos/rs_rec',data=thermodynamic_parameters['rs_rec'])
    file.create_dataset('thermos/ds_rec',data=thermodynamic_parameters['ds_rec'])
    file.create_dataset('thermos/ra_rec',data=thermodynamic_parameters['ra_rec'])
    file.create_dataset('thermos/da_rec',data=thermodynamic_parameters['da_rec'])
    file.create_dataset('thermos/rd_rec',data=thermodynamic_parameters['rd_rec'])

    file.create_dataset('thermos/rs_drag',data=cosmo.rs_drag())

    file.create_dataset('thermos/r_s',data=r_s_bg)
    file.create_dataset('thermos/r_d',data=thermos_NN["r_d"])
    file.create_dataset('thermos/e_kappa',data=thermos_NN["e_kappa"])
    file.create_dataset('thermos/conf_time_bg',data=tau_bg)
    file.create_dataset('thermos/conf_time_th',data=thermos_NN["tau"])
    file.create_dataset('thermos/g',data=thermos_NN["g"])
    file.create_dataset('thermos/g_prime',data=thermos_NN["g_prime"])
    file.create_dataset('thermos/g_reco',data=thermos_NN["g_reco"])
    file.create_dataset('thermos/g_reco_prime',data=thermos_NN["g_reco_prime"])
    file.create_dataset('thermos/g_reio',data=thermos_NN["g_reio"])
    file.create_dataset('thermos/g_reio_prime',data=thermos_NN["g_reio_prime"])

    a = np.zeros(len(tau_bg))

    for i_t,tau_i in enumerate(tau_bg):
        try:
            a[i_t] = cosmo.a_of_tau(tau_i)
        except:
            a[i_t] = a[i_t-1]

    file.create_dataset('thermos/a',data=a)

    tau_eq, a_eq, H_eq = cosmo.get_quantities_at_RM_equality()
    z_bg = cosmo.get_bg_z()
    tau_bg = cosmo.get_bg_tau()

    file.create_dataset("background/tau", data=tau_bg)
    file.create_dataset("background/z", data=z_bg)
    file.create_dataset("background/rho_b", data=rho_g)
    file.create_dataset("background/rho_g", data=rho_b)
    file.create_dataset("background/tau_eisw_lisw_50", data=cosmo.tau_of_z(50))
    file.create_dataset("background/tau_eisw_lisw_120", data=cosmo.tau_of_z(120))

    D = bg_nn["D"]
    H = bg_nn["H"]

    file.create_dataset('thermos/tau_eq',data=tau_eq)
    file.create_dataset('thermos/a_eq',data=a_eq)
    file.create_dataset('thermos/H_eq',data=H_eq)

    file.create_dataset('thermos/D',data=D)
    file.create_dataset('thermos/H',data=H)

    # keep track of the minima and maxima for normalization
    minima = {}
    maxima = {}
    for group_name, group in file.items():
        for qty_name, qty in group.items():
            # cosmological parameters are a special case:
            # they are always prefixed by 'cosmos/' in order to distinguish
            # 'cosmos/tau_reio' (i.e. the optical depth of reio) from the
            # conformal time of reio 'tau_reio'
            if group_name == "cosmos":
                key = "cosmos/" + qty_name
            else:
                key = qty_name

            minima[key] = qty[()].min()
            maxima[key] = qty[()].max()

    file.close()
    logger.debug("Wrote source function to {}".format(path))

    # prepare CLASS for next iteration
    cosmo.struct_cleanup()
    cosmo.empty()

    return minima, maxima
