import os
import multiprocessing
from tqdm import tqdm
import numpy as np
import h5py as h5
import classy
from scipy import interpolate

# At which z we want the powerspectra to be evaluated at
PK_M_AT_Z = 0.0
PK_CB_AT_Z = 2.0

def generate_spectral_data(workspace, varying_params, fixed_params, path, processes= 4, fixed_nn_only=None, use_nn=False):
    """
    Generate training/validation/testing set of size `count` by sampling the `domain` of cosmological parameters.

    The files will be placed in the `directory`.

    Optionally, the number of worker processes can be specified using `processes`;
    if not specified, the number of processes will equal the number of CPU cores.

    In addition, the global minima and maxima for each quantity over the
    entire set are written to the file "normalization.json" in `directory`
    as well.

    This function will return a Dict[str, np.ndarray] of the sampled cosmological parameters.
    """
    any_key = next(iter(varying_params))
    count = len(varying_params[any_key])

    # Prepare function arguments
    args = []
    for i in range(count):
        params = fixed_params.copy()
        cosmo_params = {key: varying_params[key][i] for key in varying_params}
        arg = (fixed_params, cosmo_params, fixed_nn_only, use_nn, workspace)
        args.append(arg)


    params = []
    cl_lenses = []
    cl_rawes = []
    for arg in tqdm(args, total=count):
       param,cl_lens, cl_raw = generate_spectra_function(arg)
       params.append(param)
       cl_lenses.append(cl_lens)
       cl_rawes.append(cl_raw)

    # TODO SG here occurs some bug
    # print("--- Start generating {} spectra with {} processes ---".format(str(len(args),str(processes))))
    # with multiprocessing.Pool(processes) as pool:
    #     for param,cl_lens, cl_raw in tqdm(pool.imap(generate_spectra_function, args), total=count):
    #         params.append(param)
    #         cl_lenses.append(cl_lens)
    #         cl_rawes.append(cl_raw)

    def create_group(f, name, data):
        g = f.create_group(name=name)
        for key in data[0].keys():
            data_column = [dictionary[key] for dictionary in data]
            if isinstance(data_column[0],str):
                data_column = [n.encode("ascii", "ignore") for n in data_column]
                g.create_dataset(key, (len(data_column),1),'S20',data = data_column)
            else:
                g.create_dataset(key, data= data_column)
    
    with h5.File(path, "w") as out:
        create_group(out, "parameter", params)
        create_group(out, "cl_lens", cl_lenses)
        create_group(out, "cl_raw", cl_rawes)


def generate_spectra_function(args):
    """
    Generate source function of cosmo_param_sets[i] and stores it into the corresponding file.
    """
    params, cosmo_params, fixed_nn_only, use_nn, workspace = args

    # create instance of the class "Class"
    cosmo = classy.Class()

    params = params.copy()
    params.update(cosmo_params)
    if use_nn==1:
        params.update({'use_nn':'yes'})
        params.update(fixed_nn_only)

    # for the purpose of testing we are interested in the linear matter power spectrum
    params['non linear'] = 'no'
    params['z_max_pk'] = '5'

    cosmo.set(**params)

    # run CLASS
    cosmo.compute()
    
    # Extract spectra
    cl_lensed = cosmo.lensed_cl()
    cl_raw    = cosmo.raw_cl()

    # Transform pk into NN k values
    k_min = cosmo.k_min()
    k_net = np.load(str(workspace.data / 'k.npy'))
    kk = np.concatenate(([k_min], k_net[(k_net>k_min)&(k_net<params["P_k_max_1/Mpc"])]))
    kk_save=np.zeros(k_net.shape, dtype=float)
    pk_save=np.zeros(k_net.shape, dtype=float)
    pk_cb_save=np.zeros(k_net.shape, dtype=float)

    pk_class, k_class = cosmo.pk_at_z(PK_M_AT_Z)
    pk_cb_class, k_class = cosmo.pk_cb_at_z(PK_CB_AT_Z)
    
    # We need to interpolate the spectra to fit the same k grid, such that we can compare them to another.
    # This is done with a cubic fit in the log space. Linear or non-log fit lead to systemtics for low k!
    class_pk_interpolation = interpolate.interp1d(np.log(k_class), np.log(pk_class), kind='cubic',fill_value='extrapolate')
    class_pk_cb_interpolation = interpolate.interp1d(np.log(k_class), np.log(pk_cb_class), kind='cubic',fill_value='extrapolate')

    for j in range(len(kk)):
        # if we did not overwrite the k-array we need to interpolate the k array
        kk_save[j]=kk[j]
        pk_save[j]=np.exp(class_pk_interpolation(np.log(kk[j])))
        pk_cb_save[j]=np.exp(class_pk_cb_interpolation(np.log(kk[j])))

    cl_lensed["pk"]=pk_save
    cl_lensed["pk_cb"]=pk_cb_save
    cl_lensed["kk"]=kk_save

    #cosmo.struct_cleanup()
    del cosmo
    return params,cl_lensed,cl_raw
