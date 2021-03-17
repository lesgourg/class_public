import os
import functools

import numpy as np
import scipy.special
import h5py as h5
import scipy.interpolate

from classynet.dependency_resolution import Graph, Evaluator

# maxima values of spherical bessel
J1_MAX_VALUE = 0.436181727
J2_MAX_VALUE = 0.306791812
INV_J1_MAX_VALUE = 1. / J1_MAX_VALUE
INV_J2_MAX_VALUE = 1. / J2_MAX_VALUE

PHI_PLUS_PSI_REFERENCE_FILE_PATH = os.path.join(
        os.path.expandvars("$CLASSNET_DATA"),
        "reference_phi_plus_psi.h5"
        )

def create_reference_spline(quantity, normalize=False):
    with h5.File(PHI_PLUS_PSI_REFERENCE_FILE_PATH, "r") as f:
        tau = f["tau"][()]
        k = f["k"][()]
        y = f[quantity][()]
        if normalize:
            y /= np.max(np.abs(y))
        # spline = scipy.interpolate.interp2d(
        #         tau,
        #         k,
        #         y,
        #         kind="cubic"
        #         )
        spline = scipy.interpolate.RectBivariateSpline(k, tau, y)
        k_eq_ref = f["k_eq"][()]
    return spline, k_eq_ref

def create_psi_minus_phi_spline(normalize=False):
    with h5.File(PHI_PLUS_PSI_REFERENCE_FILE_PATH, "r") as f:
        tau = f["tau"][()]
        k = f["k"][()]
        y = f["psi"][()] - f["phi"][()]
        if normalize:
            y /= np.max(np.abs(y))
        # spline = scipy.interpolate.interp2d(
        #         tau,
        #         k,
        #         y,
        #         kind="cubic"
        #         )
        spline = scipy.interpolate.RectBivariateSpline(k, tau, y)
    return spline

def create_reference_spline_derivative(quantity, normalize=False):
    with h5.File(PHI_PLUS_PSI_REFERENCE_FILE_PATH, "r") as f:
        tau = f["tau"][()]
        k = f["k"][()]
        # y has shape (n_k, n_tau) here; hence, for y_prime, we need
        # to take the gradient along axis 1
        y = f[quantity][()]
        y_prime = np.gradient(y, tau, axis=1)
        if normalize:
            y_prime /= np.max(np.abs(y_prime))
        # spline = scipy.interpolate.interp2d(tau, k, y_prime, kind="cubic")
        spline = scipy.interpolate.RectBivariateSpline(k, tau, y_prime)
        k_eq_ref = f["k_eq"][()]
    return spline, k_eq_ref

def create_F_spline(tau_bg, R, R_prime):
    with h5.File(PHI_PLUS_PSI_REFERENCE_FILE_PATH, "r") as f:
        tau = f["tau"][()]
        k = f["k"][()]
        phi = f["phi"][()].T
        psi = f["psi"][()].T

        phi_prime = np.gradient(phi, tau, axis=0)
        phi_prime2 = np.gradient(phi_prime, tau, axis=0)

        R = np.interp(tau, tau_bg, R)
        R_prime = np.interp(tau, tau_bg, R_prime)

        F = -phi_prime2 - (R_prime / (1 + R))[:, None] * phi_prime - k[None, :]**2 / 3 * psi

        F_spline = scipy.interpolate.interp2d(k, tau, F, kind="cubic")

        return F_spline

# phi_plus_psi_reference_spline, k_eq_ref = create_reference_spline("phi_plus_psi")
# psi_minus_phi_spline = create_psi_minus_phi_spline()
# psi_reference_spline, _ = create_reference_spline("psi", normalize=True)
# phi_reference_spline, _ = create_reference_spline("phi")
# phi_prime_reference_spline, _ = create_reference_spline_derivative("phi", normalize=True)

def step(x):
    return 0.5 + np.tanh(x / 100.) / 2.

class DummySelection:
    def __contains__(self, x):
        return True

def expand_2d(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        ret = func(*args, **kwargs)
        return np.expand_dims(ret, axis=2)
    return wrapper

class Approximations:
    graph = Graph()


    @graph.register
    def k_eq(a_eq, H_eq):
        return (a_eq * H_eq)

    @graph.register
    @expand_2d
    def sin(r_s, k):
        return np.sin(np.outer(r_s, k))

    @graph.register
    @expand_2d
    def cos(r_s, k):
        return np.cos(np.outer(r_s, k))

    @graph.register
    def j1_j2(k, tau, tau_rec):
        return libbessel.compute_bessels(k, tau - tau_rec)

    @graph.register
    @expand_2d
    def j1(j1_j2):
        j1, _ = j1_j2
        return j1

    @graph.register
    @expand_2d
    def j2(j1_j2):
        _, j2 = j1_j2
        return j2

    @graph.register
    @expand_2d
    def damping(k, k_d):
        return np.exp(-(k[None, :] / k_d[:, None])**2)

    @graph.register
    @expand_2d
    def t0_reco_approx_1(g_reco_prime, sin, damping, k):
        return g_reco_prime[:, None] * sin[..., 0] * damping[..., 0] / k[None, :]

    @graph.register
    @expand_2d
    def t0_reco_approx_2(g_reco, cos, damping):
        return g_reco[:, None] * cos[..., 0] * damping[..., 0]

    @graph.register
    @expand_2d
    def t0_reco_approx_3(g_reco_prime, cos, damping, k):
        return g_reco_prime[:, None] * cos[..., 0] * damping[..., 0] / k[None, :]

    @graph.register
    @expand_2d
    def t0_reco_approx_4(g_reco, sin, damping):
        return g_reco[:, None] * sin[..., 0] * damping[..., 0]

    @graph.register
    @expand_2d
    def t0_reco_approx_5(g_reco, tau, k, k_eq):
        return g_reco[:, None] * psi_reference_spline(k * k_eq_ref / k_eq[0], tau).T

    @graph.register
    @expand_2d
    def t0_reco_approx_6(g_reco_prime, cos, damping):
        return g_reco_prime[:, None] * cos[..., 0] * damping[..., 0]

    @graph.register
    @expand_2d
    def t0_reco_approx_7(g_reco, sin, k, damping):
        return g_reco[:, None] * sin[..., 0] / k[None, :] * damping[..., 0]

    @graph.register
    @expand_2d
    def t0_reco_approx_8(g_reco, cos, k, damping):
        return g_reco[:, None] * cos[..., 0] / k[None, :] * damping[..., 0]

    @graph.register
    def si_ci(k, r_s):
        si, ci = scipy.special.sici(np.outer(r_s, k))
        si -= np.pi / 2
        return si, ci

    @graph.register
    def si(si_ci):
        si, _ = si_ci
        return si

    @graph.register
    def ci(si_ci):
        _, ci = si_ci
        return ci

    @graph.register
    def t0_reco_basis(tau, tau_rec, k, r_s, sin, cos, damping, si, ci):
        # raise RuntimeError("Should not be called")
        get_basis_sin = lambda p: k[None, :]**p * sin[..., 0] * damping[..., 0] / k[None, :]
        get_basis_cos = lambda p: k[None, :]**p * cos[..., 0] * damping[..., 0]
        max_power = 4
        powers = np.arange(max_power + 1)
        basis_sin = [get_basis_sin(p) for p in powers]
        basis_cos = [get_basis_cos(p) for p in powers]

        sid = si * damping[..., 0]
        cid = ci * damping[..., 0]

        r_s = r_s[:, None]
        x = k[None, :] * r_s

        xp = np.power(x[None, ...], powers[:, None, None])
        rsp = np.power(r_s[None, ...], powers[:, None, None])

        irs = 1 / r_s

        basis_si = [
                sid,
                x * sid + basis_cos[0] + basis_sin[0] * irs,
                xp[2] * sid + r_s * basis_cos[1] + basis_sin[1],
                xp[3] * sid + rsp[2] * basis_cos[2] - 2 * basis_cos[0] + r_s * basis_sin[2] - 6 * basis_sin[0] * irs,
                xp[4] * sid + rsp[3] * basis_cos[3] - 2 * r_s * basis_cos[1] + rsp[2]  * basis_sin[3] - 6 * basis_sin[1],
                ]

        basis_ci = [
                cid - basis_sin[0] / r_s,
                x * cid - basis_sin[1],
                xp[2] * cid - r_s * basis_sin[2] + 2 * basis_sin[0] * irs + basis_cos[0],
                xp[3] * cid - rsp[2] * basis_sin[3] + 2 * basis_sin[1] + r_s * basis_cos[1],
                xp[4] * cid - rsp[3] * basis_sin[4] + 2 * r_s * basis_sin[2] + rsp[2] * basis_cos[2] - 24 * basis_sin[0] * irs - 6 * basis_cos[0],
                ]

        # We do not want the normalization to depend on what tau sampling is used for
        # network training and evaluation. Hence, we should always normalize w.r.t. to
        # the basis functions' maxima at reco.

        # Make sure that reco is actually contained in the tau array, otherwise, this is pointless
        tau_rec = tau_rec[0]
        assert tau.min() < tau_rec
        assert tau.max() > tau_rec

        basis = np.dstack(basis_sin + basis_cos + basis_si + basis_ci)
        # Create a spline to get the basis functions at rec precisely (so that there are no rounding artifacts, etc.)
        spline = scipy.interpolate.interp1d(tau, basis, axis=0)
        basis_at_rec = spline(tau_rec)
        # basis_at_rec = basis[np.argmin(np.abs(tau - tau_rec))]

        basis_normalized = basis / np.abs(basis_at_rec).max(axis=0)

        return basis_normalized

    @graph.register
    def t0_reco_basis_with_g(g_reco, t0_reco_basis):
        # raise RuntimeError("Should not be called")
        return g_reco[:, None, None] * t0_reco_basis

    @graph.register
    def t0_reco_basis_with_g_prime(g_reco_prime, t0_reco_basis):
        # raise RuntimeError("Should not be called")
        return g_reco_prime[:, None, None] * t0_reco_basis

    @graph.register
    @expand_2d
    def t2_reco_approx_1(g_reco, sin, damping):
        return g_reco[:, None] * sin[..., 0] * damping[..., 0]

    @graph.register
    @expand_2d
    def t2_reco_approx_1(g_reco, sin, damping, k):
        return g_reco[:, None] * sin[..., 0] * k[None, :]**2 * damping[..., 0]

    @graph.register
    def interp_param(tau, tau_rec):
        return np.expand_dims(step(tau - 4 * tau_rec), axis=1)

    @graph.register
    @expand_2d
    def cosj1(cos, j1, interp_param):
        t = interp_param
        return (1 - t) * cos[..., 0] + t * j1[..., 0]

    @graph.register
    @expand_2d
    def sinj2(sin, j2, interp_param):
        t = interp_param
        return (1 - t) * sin[..., 0] + t * j2[..., 0]

    @graph.register
    @expand_2d
    def reference(tau, k, k_eq):
        return phi_plus_psi_reference_spline(k * k_eq_ref / k_eq[0], tau).T

    @graph.register
    def D_prime(D, tau):
        return np.gradient(D, tau)

    @graph.register
    @expand_2d
    def t0_reio_approx_1(g_reio_prime, D_prime, reference):
        return g_reio_prime[:, None] * D_prime[:, None] * reference[:, :, 0]

    @graph.register
    def phi_prime(tau, k, k_eq):
        return phi_prime_reference_spline(k * k_eq_ref / k_eq[0], tau).T

    @graph.register
    @expand_2d
    def psi(tau, k, k_eq):
        return psi_reference_spline(k * k_eq_ref / k_eq[0], tau).T

    @graph.register
    def psi_minus_phi(tau, k, k_eq):
        return psi_minus_phi_spline(k * k_eq_ref / k_eq[0], tau).T

    @graph.register
    def raw_g_reco(g_reco):
        return g_reco

    @graph.register
    def raw_tau(tau):
        return tau

    @graph.register
    def raw_g_reco_prime(g_reco_prime):
        return g_reco_prime

    @graph.register
    def R(rho_b, rho_g):
        return 3 * rho_b / (4 * rho_g)

    @graph.register
    def raw_tau(tau):
        return tau

    evaluator = graph.evaluator()

def add_approximations(inputs, k, tau, selection=None):

    initial_cache = inputs.copy()

    # approximations etc. expect k vector, not expanded in tau axis
    initial_cache["k"] = k

    if selection is not None:
        ret = Approximations.evaluator.evaluate_many(selection, cache=initial_cache)
        # remove the vector k (so the expanded k in `inputs` doesn't get overwritten)
        if "k" in ret:
            del ret["k"]
    else:
        ret = Approximations.evaluator.evaluate_all(cache=initial_cache)

    inputs.update(ret)
