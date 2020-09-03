import glob
import os
import numpy as np
from tqdm import tqdm

from parameter_sampling import EllipsoidDomain

COL = {
    "omega_b": 2,
    "omega_cdm": 3,
    "tau_reio": 5,
    "H0": 29
}

# path = "/scratch/work/samaras/planck_stuff/COM_CosmoParams_fullGrid_R3.01/base/plikHM_TTTEEE_lowl_lowE_lensing/"
# filenames = ["base_plikHM_TTTEEE_lowl_lowE_lensing_{}.txt".format(i) for i in [1, 2, 3, 4]]

# path = "/scratch/work/samaras/planck_stuff/COM_CosmoParams_fullGrid_R3.01/base/plikHM_TTTEEE_lowE/"
# filenames = ["base_plikHM_TTTEEE_lowE_{}.txt".format(i) for i in [1, 2, 3, 4]]

path = "/scratch/work/samaras/planck_stuff/COM_CosmoParams_fullGrid_R3.01/base/plikHM_TTTEEE_lowl_lowE/"
filenames = ["base_plikHM_TTTEEE_lowl_lowE_{}.txt".format(i) for i in [1, 2, 3, 4]]

data = []

print("loading planck chain...")
for fname in filenames:
    p = os.path.join(path, fname)
    data.append(np.genfromtxt(p))
data = np.concatenate(data)
print("done.")

# omega_b = data[:, COL["omega_b"]]
# omega_cdm = data[:, COL["omega_cdm"]]
# tau_reio = data[:, COL["tau_reio"]]
# H0 = data[:, COL["H0"]]
# print(np.mean(H0))
# bf_index = data[:, 1].argmax()
# print("H0 bf:", H0[bf_index])

print("creating param_dicts...")
param_dicts = [{name: row[COL[name]] for name in COL} for row in data]
print("done.")

### Planck + BAO + Pantheon parameter ellipsoid ###

pnames = [
    'omega_b',
    'omega_cdm',
    'h',
    'tau_reio',
    'w0_fld',
    'wa_fld',
    'N_ur',
    'omega_ncdm',
    'Omega_k']

pnames_ext = [
    'w0_fld',
    'wa_fld',
    'N_ur',
    'omega_ncdm',
    'Omega_k']

DATA_DIR = "/home/samaras/CLASSnet_HPC/data"

domain = EllipsoidDomain(
    bestfit_path=os.path.join(DATA_DIR, "lcdm_11p_sn.bestfit"),
    covmat_path=os.path.join(DATA_DIR, "lcdm_11p_sn.covmat"),
    pnames=pnames,
    sigma=7
)

# we need to augment param_dicts to contain the full 9 parameters
# we do this by setting all extended parameters to their bestfit value
assert len(domain.best_fit) == 9

extra = {name: domain.best_fit[domain.index(name)] for name in pnames_ext}

print("checking if inside ellipsoid...")
inside = 0
for row in tqdm(param_dicts):
    pdict = row.copy()
    pdict.update(extra)
    H0 = pdict.pop("H0")
    pdict["h"] = H0 / 100
    if domain.contains(pdict):
        inside += 1
print("done.")

inside_fraction = inside / len(param_dicts)
print("Fraction of Planck baseline points inside (Planck+BAO+Pantheon) ellipsoid:", inside_fraction)
