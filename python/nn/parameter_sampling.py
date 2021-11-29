import os
import numpy as np
import scipy.stats as ss
from classynet.lhs import lhs, lhs_float
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

class EllipsoidDomain(ParamDomain):

    @staticmethod
    def from_paths(bestfit_path, covmat_path, pnames, sigma_train, sigma_validate):
        _bestfit_names, best_fit = load_montepython_bestfit(bestfit_path, pnames)
        assert _bestfit_names == pnames
        covmat, inv_covmat = load_montepython_covmat(covmat_path, pnames)
        return EllipsoidDomain(
            best_fit=best_fit,
            covmat=covmat,
            inv_covmat=inv_covmat,
            pnames=pnames,
            sigma_train=sigma_train,
            sigma_validate=sigma_validate,
        )

    def __init__(self, best_fit, covmat, inv_covmat, pnames, sigma_train, sigma_validate=None):
        self.best_fit = best_fit
        self.covmat = covmat
        self.inv_covmat = inv_covmat
        self.pnames = pnames
        self.sigma_train = sigma_train
        self.sigma_validate = sigma_validate if sigma_validate is not None else self.sigma_train

        self.ndf = len(self.pnames)

    def index(self, name):
        return self.pnames.index(name)

    def sample(self, count, sigma, tol=50, maxiter=100):
        fraction = 1.0
        for _ in range(maxiter):
            if fraction > 0:
                new_count = int(count / fraction)
            else:
                new_count = 5 * new_count
            samples = self._sample(new_count, sigma=sigma)
            print("-----------------------------------------------")
            if abs(len(samples) - count) < tol:
                return samples
            fraction = len(samples) / new_count
        raise ValueError("Couldn't get {} samples with tol={}".format(count, tol))

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
        print("samples.shape", samples.shape)

        DISABLE_OMEGA_K = False
        if DISABLE_OMEGA_K:
            # !!!!! TODO DANGER !!!!!!
            # MAKE SURE THIS IS NOT ENABLED!
            # !!!!! TODO DANGER !!!!!!
            samples[:, self.index("Omega_k")] = 0.0
            self.best_fit[self.index("Omega_k")] = 0.0
            ndf = self.ndf - 1
        else:
            ndf = self.ndf

        inside_ellipsoid = is_inside_ellipsoid(
            self.inv_covmat,
            samples, self.best_fit,
            ndf, sigma)
        ratio_inside = inside_ellipsoid.sum() / len(samples)
        print("fraction of points inside ellipsoid:", ratio_inside)

        tau_large_enough = samples[:, self.index("tau_reio")] > 0.004
        count_tau = tau_large_enough.sum()
        ratio_tau = count_tau / len(samples)
        print("fraction of kept points (tau_reio only):", ratio_tau)

        N_ur_positive = samples[:, self.index("N_ur")] > 0
        ratio_N_ur = N_ur_positive.sum() / len(samples)
        print("fraction of kept points (N_ur only):", ratio_N_ur)

        omega_ncdm_large_enough = samples[:, self.index("omega_ncdm")] > OMEGA_NCDM_MIN
        ratio_ncdm = omega_ncdm_large_enough.sum() / len(samples)
        print("fraction of kept points (omega_ncdm only):", ratio_ncdm)

        fld_consistent = samples[:, self.index("w0_fld")] + samples[:, self.index("wa_fld")] < -1./3.
        count_fld = fld_consistent.sum()
        ratio_fld = count_fld / len(samples)
        print("fraction of kept points (fld only):", ratio_fld)

        # wa_bound = samples[:, self.index("wa_fld")] <= 0.0
        # count_wa = wa_bound.sum()
        # ratio_wa = count_wa / len(samples)
        # print("fraction of kept points (wa < 0 only):", ratio_wa)
        # inside_mask = inside_ellipsoid & tau_large_enough & fld_consistent & omega_ncdm_large_enough & wa_bound

        inside_mask = inside_ellipsoid \
            & tau_large_enough \
            & fld_consistent \
            & omega_ncdm_large_enough \
            & N_ur_positive
        count_inside = inside_mask.sum()
        ratio_inside = count_inside / len(samples)
        print("count_inside:", count_inside)
        print("fraction of kept points:", ratio_inside)

        samples = samples[inside_mask]

        return samples

    def sample_save(self, training_count, validation_count, path):
        def create_group(f, name, count, sigma):
            samples = self.sample(count, sigma=sigma)
            print("Saving group '{}' of {} samples to {}".format(name, len(samples), path))
            g = f.create_group(name=name)
            for name, column in zip(self.pnames, samples.T):
                g.create_dataset(name, data=column)

        with h5.File(path, "w") as out:
            create_group(out, "training", training_count, sigma=self.sigma_train)
            create_group(out, "validation", validation_count, sigma=self.sigma_validate)

    def parameter_names(self):
        return list(self.pnames)

    def contains(self, parameters, validate=True):
        """
        parameters is a dict of CLASS parameters.
        NOTE: this assumes that all network parameters (i.e. all of self.pnames)
        are present in parameters.
        """
        if parameters["tau_reio"] <= 0.004:
            return False, 1001
        if parameters["w0_fld"] + parameters["wa_fld"] > 0:
            return False, 1002
        if parameters["omega_ncdm"] <= OMEGA_NCDM_MIN:
            return False, 1003
        if parameters["N_ur"]<0:
            return False, 1004
        cosmo_params = np.array([parameters[name] for name in self.pnames])[None, :]
        sigma = self.sigma_validate if validate else self.sigma_train
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
            "sigma_validate": self.sigma_validate,
        }
        with open(path, "w") as out:
            json.dump(d, out)

    @staticmethod
    def load(path):
        with open(path, "r") as src:
            data = json.load(src)
        return EllipsoidDomain(
            best_fit       = data["best_fit"],
            covmat         = data["covmat"],
            inv_covmat     = data["inv_covmat"],
            pnames         = data["pnames"],
            sigma_train    = data["sigma_train"],
            sigma_validate = data["sigma_validate"],
        )

class DefaultParamDomain(ParamDomain):
    def __init__(self, covmat_planck_path, sigma):
        self.sigma = sigma
        self.planck_names = ["omega_b", "omega_cdm", "tau_reio", "H0"]
        self.ndf = len(self.planck_names)
        self.best_fit_planck = np.array([
            # omega_b
            0.02236,
            # omega_cdm
            0.1202,
            # tau_reio
            0.0544,
            # H0
            67.27,
        ])
        self.other_names = ["w0_fld", "wa_fld", "N_ur", "omega_ncdm", "Omega_k"]
        self.other_bounds = np.array([
            # w0_fld
            [-1.5, -0.5],
            # wa_fld
            [-5.0, -0.4],
            # N_ur
            [ 0.0,  0.1],
            # omega_ncdm:
            # lower bound is minimum value that doesn't crash class
            # upper bound is such that m_tot ~ 0.9eV
            [8.147986e-5, 0.00966],
            # Omega_k
            [-0.03, 0.03],
        ])
        self.names = self.planck_names + self.other_names
        self.ndf_total = len(self.names)
        self.covmat, self.inv_covmat = load_montepython_covmat(covmat_planck_path, self.planck_names)

    def index(self, name):
        return self.names.index(name)

    def sample(self, count, tol=100, maxiter=10):
        fraction = 1.0
        for _ in range(maxiter):
            new_count = int(count / fraction)
            samples = self._sample(new_count)
            print("-----------------------------------------------")
            if abs(len(samples) - count) < tol:
                return samples
            fraction = len(samples) / new_count
        raise ValueError("Couldn't get {} samples with tol={}".format(count, tol))


    def _sample(self, count):
        deltachi2 = get_delta_chi2(self.ndf, self.sigma)
        bbox_planck = find_bounding_box(self.best_fit_planck, self.covmat, deltachi2)

        # combined bounding box
        bbox = np.concatenate((bbox_planck, self.other_bounds))
        lower, upper = bbox[:, 0], bbox[:, 1]

        from time import perf_counter
        start = perf_counter()
        # samples = lhs(self.ndf_total, samples=count)
        samples = lhs_float(self.ndf_total, count).T
        elapsed = perf_counter() - start

        print("lhs took {}s".format(elapsed))
        samples = lower[None, :] + samples * (upper - lower)[None, :]
        print("samples.shape", samples.shape)
        samples_planck = samples[:, :len(self.planck_names)]


        inside_planck_ellipsoid = is_inside_ellipsoid(
            self.inv_covmat,
            samples_planck, self.best_fit_planck,
            self.ndf, self.sigma)
        count_planck = inside_planck_ellipsoid.sum()
        ratio_planck = count_planck / len(samples)
        print("fraction of kept points (planck only):", ratio_planck)

        tau_large_enough = samples_planck[:, self.index("tau_reio")] > 0.004
        count_tau = tau_large_enough.sum()
        ratio_tau = count_tau / len(samples)
        print("fraction of kept points (tau_reio only):", ratio_tau)

        fld_consistent = samples[:, self.index("w0_fld")] + samples[:, self.index("wa_fld")] < -1./3.
        count_fld = fld_consistent.sum()
        ratio_fld = count_fld / len(samples)
        print("fraction of kept points (fld only):", ratio_fld)

        inside_mask = inside_planck_ellipsoid & tau_large_enough & fld_consistent
        count_inside = inside_mask.sum()
        ratio_inside = count_inside / len(samples)
        print("count_inside:", count_inside)
        print("fraction of kept points:", ratio_inside)

        samples = samples[inside_mask]

        return samples

    def sample_save(self, training_count, validation_count, path):
        def create_group(f, name, count):
            samples = self.sample(count)
            print("Saving group '{}' of {} samples to {}".format(name, len(samples), path))
            g = f.create_group(name=name)
            for name, column in zip(self.names, samples.T):
                g.create_dataset(name, data=column)

        with h5.File(path, "w") as out:
            create_group(out, "training", training_count)
            create_group(out, "validation", validation_count)

    def contains(self, parameters):
        raise NotImplementedError

    def save(self, path):
        d = {
            "fields": list(self.names),
            "best_fit_planck": list(self.best_fit_planck),
            "covmat": [list(row) for row in self.covmat],
            "inv_covmat": [list(row) for row in self.inv_covmat],
            "other_bounds": list(self.other_bounds)
        }
        with open(path, "w") as out:
            json.dump(d, out)

    def load(self, path):
        raise NotImplementedError


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
    delta_chi2 = get_delta_chi2(df, sigma)
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
