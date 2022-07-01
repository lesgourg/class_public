import torch


def get_fields(x, fields):
    return torch.stack([x[s] for s in fields], 1)

cosmo_names = ["omega_b", "omega_cdm", "h", "tau_reio", "w0_fld", "wa_fld", "N_ur", "omega_ncdm", "Omega_k"]
INPUTS_COSMO = ["cosmos/" + name for name in cosmo_names]

def get_inputs_cosmo(x):
    return get_fields(x, INPUTS_COSMO)

def get_inputs_tau_reco(x):
    fields = ("tau_relative_to_reco", "g_reco", "g_reco_prime", "e_kappa")
    return get_fields(x, fields)


def get_inputs_tau_reio(x):
    fields = ("tau_relative_to_reio", "g_reio", "g_reio_prime", "e_kappa")
    return get_fields(x, fields)

def get_inputs_tau_isw(x):
    return get_fields(x, ("tau", "e_kappa"))

########################## LOSS FUNCTIONS ####################################


def mse_truncate_(k, k_min):
    import numpy as np
    # index of first member in `k` which is larger than k_min
    k_min_idx = np.searchsorted(k.cpu().numpy(), k_min.cpu().item())
    # shift by one to get the index of the first member in `k` which
    # is smaller than `k_min` (and make sure we don't fall out of the array)
    k_min_idx = max(k_min_idx - 1, 0)
    def loss(prediction, truth):
        # only keep points starting from k_min_idx
        return torch.mean((prediction - truth)[:, k_min_idx:, ...]**2)
    return loss

def mse_():
    def loss(prediction, truth):
        return torch.mean((prediction - truth)**2)
    return loss

def mse_truncate(k, k_min):
    return mse_truncate_(k, k_min)

def mse_rel_():
    def loss(prediction, truth):
        return torch.mean(((prediction - truth) / truth)**2)
    return loss

def mse_rel_truncate_(k, k_min):
    mask = k >= k_min
    def loss(prediction, truth):
        return torch.mean(((prediction - truth) / truth)[:, mask, ...]**2)
    return loss

def mse_rel_truncate(k, k_min):
    return mse_rel_truncate_(k, k_min)
