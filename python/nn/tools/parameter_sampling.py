from genericpath import exists
import os
import time
import numpy as np
import scipy.stats as ss
from classynet.tools.lhs import lhs, lhs_float
from tabulate import tabulate
import json
import h5py as h5

class ParamDomain:
    def sample(self, count):
        raise NotImplementedError

    def contains(self, parameters):
        raise NotImplementedError

    def save(self, path):
        raise NotImplementedError

    def save_h5(self, path):
        raise NotImplementedError

    def load(self, path):
        raise NotImplementedError

#OMEGA_NCDM_MIN = 8.147986e-5
OMEGA_NCDM_MIN = 1.70698158e-5
W0WA_BOUND = -1./3.

class EllipsoidDomain(ParamDomain):

    @staticmethod
    def from_paths(workspace,
        pnames, 
        bestfit_path, 
        covmat_path, 
        sigma_train = 6, 
        sigma_validation = 5,
        sigma_test = 5,
        ):
        _bestfit_names, best_fit = load_montepython_bestfit(bestfit_path, pnames)
        assert _bestfit_names == pnames
        covmat, inv_covmat = load_montepython_covmat(covmat_path, pnames)
        return EllipsoidDomain(
            workspace,
            best_fit=best_fit,
            covmat=covmat,
            inv_covmat=inv_covmat,
            pnames=pnames,
            sigma_train=sigma_train,
            sigma_validation=sigma_validation,
            sigma_test=sigma_test,
        )

    def __init__(self, workspace, best_fit, covmat, inv_covmat, pnames, sigma_train=6, sigma_validation=5, sigma_test=5):
        self.workspace = workspace
        self.best_fit = best_fit
        self.covmat = covmat
        self.inv_covmat = inv_covmat
        self.pnames = pnames
        self.sigma_train = sigma_train
        self.sigma_validation = sigma_validation 
        self.sigma_test = sigma_test 

        self.ndf = len(self.pnames)

    def index(self, name):
        return self.pnames.index(name)

    def sample(self, count, sigma, tol=0.1, maxiter=100):
        fraction = 1.0
        for _ in range(maxiter):
            if fraction > 0:
                new_count = int(count / fraction)
            else:
                new_count = 5 * new_count
            samples = self._sample(new_count, sigma=sigma)
            print("-----------------------------------------------")
            if abs(len(samples) - count)/count < tol:
                #wait a second to ensure the get a new seed is created for following dataset
                time.sleep(1)
                return samples
            fraction = len(samples) / new_count
        raise ValueError("Couldn't get {} samples with relative tol={}".format(count, tol))

    def _sample(self, count, sigma):
        deltachi2 = get_delta_chi2(self.ndf, sigma)
        bbox = find_bounding_box(self.best_fit, self.covmat, deltachi2)
        # lower and upper boundaries of bounding box define lhs domain
        lower, upper = bbox[:, 0], bbox[:, 1]

        from time import perf_counter
        start = perf_counter()
        # samples = lhs(self.ndf, samples=count)
        samples = lhs_float(self.ndf, count).T
        assert samples.shape == (count, self.ndf)
        elapsed = perf_counter() - start

        print("lhs took {}s".format(elapsed))
        samples = lower[None, :] + samples * (upper - lower)[None, :]

        inside_ellipsoid = is_inside_ellipsoid(
            self.inv_covmat,
            samples, self.best_fit,
            self.ndf, sigma)

        ratio_inside = inside_ellipsoid.sum() / len(samples)
        print("fraction of points inside ellipsoid:", ratio_inside)

        # eff_mask will allow 
        eff_mask = np.array(inside_ellipsoid,dtype=np.bool)

        if "tau_reio" in self.pnames:
            tau_large_enough = samples[:, self.index("tau_reio")] > 0.004
            count_tau = tau_large_enough.sum()
            ratio_tau = count_tau / len(samples)
            eff_mask = eff_mask & tau_large_enough
            print("fraction of kept points (tau_reio only):", ratio_tau)
        
        if "N_ur" in self.pnames:
            N_ur_positive = samples[:, self.index("N_ur")] > 0
            ratio_N_ur = N_ur_positive.sum() / len(samples)
            eff_mask = eff_mask & N_ur_positive
            print("fraction of kept points (N_ur only):", ratio_N_ur)

        if "omega_ncdm" in self.pnames:
            omega_ncdm_large_enough = samples[:, self.index("omega_ncdm")] > OMEGA_NCDM_MIN
            ratio_ncdm = omega_ncdm_large_enough.sum() / len(samples)
            eff_mask = eff_mask & omega_ncdm_large_enough
            print("fraction of kept points (omega_ncdm only):", ratio_ncdm)

        if "w0_fld" in self.pnames:
            fld_consistent = samples[:, self.index("w0_fld")] + samples[:, self.index("wa_fld")] < W0WA_BOUND 
            count_fld = fld_consistent.sum()
            ratio_fld = count_fld / len(samples)
            eff_mask = eff_mask & fld_consistent
            print("fraction of kept points (fld only):", ratio_fld)


        count_inside = eff_mask.sum()
        ratio_inside = count_inside / len(samples)
        print("count_inside:", count_inside)
        print("fraction of kept points:", ratio_inside)

        samples = samples[eff_mask]

        return samples

    def sample_save(self, training_count=0, validation_count=0, test_count=0, file_name = 'parameter_sample'):
        def create_group(f, name, count, sigma):
            samples = self.sample(count, sigma=sigma)
            print("Saving group '{}' of {}".format(name, len(samples)))
            g = f.create_group(name=name)
            for name, column in zip(self.pnames, samples.T):
                g.create_dataset(name, data=column)

        # Write and sample
        if training_count!=0:
            os.makedirs(self.workspace.training_data, exist_ok=True)
            with h5.File(self.workspace.training_data / '{}.h5'.format(file_name), "w") as out:
                create_group(out, "training", training_count, sigma=self.sigma_train)
        if validation_count!=0:
            os.makedirs(self.workspace.validation_data, exist_ok=True)
            with h5.File(self.workspace.validation_data / '{}.h5'.format(file_name), "w") as out:
                create_group(out, "validation", validation_count, sigma=self.sigma_validation)
        if test_count!=0:
            os.makedirs(self.workspace.test_data, exist_ok=True)
            with h5.File(self.workspace.test_data / '{}.h5'.format(file_name), "w") as out:
                create_group(out, "test", test_count, sigma=self.sigma_test)

    def parameter_names(self):
        return list(self.pnames)

    def contains(self, parameters, validate=True):
        """
        parameters is a dict of CLASS parameters.
        NOTE: this assumes that all network parameters (i.e. all of self.pnames)
        SG: TODO CHANGE THIS!
        are present in parameters.
        """
        if "tau_reio" in parameters.keys():
            if parameters["tau_reio"] <= 0.004:
                return False, 1001
        if "w0_fld" in parameters.keys():
            if parameters["w0_fld"] + parameters["wa_fld"] > W0WA_BOUND:
                return False, 1002
        else:
            parameters["w0_fld"]=-1
        if "wa_fld" in parameters.keys():
            pass
        else:
            parameters["wa_fld"]=0
        if "Omega_k" in parameters.keys():
            pass
        else:
            parameters["Omega_k"]=0
        if "omega_ncdm" in parameters.keys():
            if parameters["omega_ncdm"] <= OMEGA_NCDM_MIN:
                return False, 1003
        else:
            parameters["omega_ncdm"]=0.06/93.1
        if "N_ur" in parameters.keys():
            if parameters["N_ur"]<0:
                return False, 1004
        else:
            parameters["N_ur"]=0.00641
        cosmo_params = np.array([parameters[name] for name in self.pnames])[None, :]
        sigma = self.sigma_validation if validate else self.sigma_train
        inside, delta_chi2 = is_inside_ellipsoid(self.inv_covmat, cosmo_params,
                                     self.best_fit, self.ndf, sigma, return_delta_chi2=True)

        return inside[0], delta_chi2[0]

    def save(self, path):
        d = {
            "best_fit":       list(self.best_fit),
            "covmat":         [list(row) for row in self.covmat],
            "inv_covmat":     [list(row) for row in self.inv_covmat],
            "pnames":         list(self.pnames),
            "sigma_train":    self.sigma_train,
            "sigma_validation": self.sigma_validation,
            "sigma_test": self.sigma_test,
        }
        with open(path, "w") as out:
            json.dump(d, out)

    @staticmethod
    def load(workspace, path):
        with open(path, "r") as src:
            data = json.load(src)
        return EllipsoidDomain(
            workspace,
            best_fit       = data["best_fit"],
            covmat         = data["covmat"],
            inv_covmat     = data["inv_covmat"],
            pnames         = data["pnames"],
            sigma_train    = data["sigma_train"],
            sigma_validation = data["sigma_validation"],
            sigma_test = data["sigma_test"],
        )

def load_montepython_file(fname):
    with open(fname) as f:
        line = f.readline().strip()
        # remove leading '#'
        line = line[1:].strip()
        names = [name.strip() for name in line.split(",")]
    return names, np.genfromtxt(fname, skip_header=1)

def load_montepython_bestfit(fname, parameters=None):
    """
    Load bestfit from given `fname` and return as dict of `name: value`.
    If `parameters` is not None, it must be a list of parameter names
    to be filtered in the result.
    """
    names, data = load_montepython_file(fname)
    if not parameters:
        return names, data
    else:
        d = dict(zip(names, data))
        return parameters, np.array([d[k] for k in parameters])

def load_montepython_covmat(fname, parameters=None):
    """
    load covariance matrix from `fname` and return (covmat, inv_covmat).
    If `bool(parameters)`, take the submatrices defined by the list of names
    `parameters`.
    """
    names, covmat = load_montepython_file(fname)
    inv_covmat = np.linalg.inv(covmat)

    if not parameters:
        return covmat, inv_covmat

    indices = [names.index(key) for key in parameters]

    submask = np.ix_(indices, indices)
    covmat_reduced = covmat[submask]
    inv_covmat_reduced = np.linalg.inv(covmat_reduced)

    return covmat_reduced, inv_covmat_reduced

def get_delta_chi2(df, sigma):
    prob = 1 - 2 * ss.norm.sf(sigma)
    return ss.chi2.ppf(prob,df)

def is_inside_ellipsoid(inv_covmat, params, bestfit, df, sigma, return_delta_chi2 = False):
    """
    params: (n, k)
    bestfit (k,)
    df: float
    sigma: float
    """
    #return([True],[0])
    delta_chi2 = get_delta_chi2(df, sigma)#*1000
    d = params - bestfit
    if return_delta_chi2 == False:
        return (d.dot(inv_covmat) * d).sum(axis=1) < delta_chi2
    else:
        return (d.dot(inv_covmat) * d).sum(axis=1) < delta_chi2, (d.dot(inv_covmat) * d).sum(axis=1)


def is_inside_ellipsoid_(inv_covmat, params, bestfit, mask, df, sigma):
    """
    params: (n, k)
    bestfit (k,)
    df: float
    sigma: float
    """
    delta_chi2 = get_delta_chi2(df, sigma)
    d = params - bestfit
    d[~mask] = 0.0
    return (d.dot(inv_covmat) * d).sum(axis=1) < delta_chi2

def find_bounding_box(bestfit, covmat, chi2):
    """
    return a (len(bestfit), 2) array of lower and upper
    bounds for the bounding box of the ellipsoid
    defined by covmat
    """
    d = np.diag(covmat)
    offsets = np.sqrt(chi2 * d)
    return np.stack([bestfit - offsets, bestfit + offsets], axis=1)
