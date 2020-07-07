import numpy as np
from .containers import InputContainer, TargetContainer

class SourceFileAdapter:
    def __init__(self, source_file):
        self.source_file = source_file

    def get_inputs(self):
        f = self.source_file

        container = InputContainer()

        thermos = f["thermos"]
        background = f["background"]

        # tau sampling of the background quantities
        container.tau_bg = thermos["conf_time_bg"][()]
        # tau sampling of the thermodynamics quantities
        container.tau_th = thermos["conf_time_th"][()]

        # Load some characteristic times
        container.scalars.tau_rec = thermos["tau_rec"][()]
        container.scalars.tau_reio = thermos["tau_reio"][()]
        container.scalars.tau_eq = thermos["tau_eq"][()]

        container.scalars.tau_eisw_lisw_50 = background["tau_eisw_lisw_50"]
        container.scalars.tau_eisw_lisw_120 = background["tau_eisw_lisw_120"]

        container.scalars.rs_drag = thermos["rs_drag"][()]

        # values at equality
        container.scalars.a_eq = thermos['a_eq'][()]
        container.scalars.H_eq = thermos['H_eq'][()]

        for qty in ("rs_rec", "ds_rec", "ra_rec", "da_rec", "rd_rec"):
            value = thermos[qty][()]
            container.scalars[qty] = value

        container.background.r_s = thermos["r_s"][()]
        container.background.D = thermos["D"][()]

        container.background.rho_b = background["rho_b"]
        container.background.rho_g = background["rho_g"]

        container.thermo.e_kappa = thermos["e_kappa"][()]
        container.thermo.r_d = thermos["r_d"][()]

        container.thermo.g = thermos["g"][()]
        container.thermo.g_prime = thermos["g_prime"][()]
        container.thermo.g_reco = thermos["g_reco"][()]
        container.thermo.g_reco_prime = thermos["g_reco_prime"][()]
        container.thermo.g_reio = thermos["g_reio"][()]
        container.thermo.g_reio_prime = thermos["g_reio_prime"][()]


        container.tau_source = f["sampling/tau"][()]

        container.scalars.k_min = f["sampling/k"][()][0]

        ###################################################################
        # 2) process cosmological_parameters
        ###################################################################
        for key, value in f["cosmos"].items():
            # Prefix cosmological parameters with 'cosmos/' to distinguish
            # the cosmological parameter "tau_rec" (i.e. the optical depth
            # to reionization) from the thermodynamical parameter "tau_rec"
            # (the conformal time of reionization)
            # This is done in order to be consistent with the naming scheme
            # of input parameters of CLASS
            container.cosmos["raw_cosmos/" + key] = np.array(value)
            container.cosmos["cosmos/" + key] = np.array(value)

        return container

    def get_outputs(self, selection=None):
        container = TargetContainer()

        f = self.source_file

        container.k_source = f["sampling/k"][()]
        container.tau_source = f["sampling/tau"][()]

        group = f["sources"]

        # Extract the source function subset defined by
        # the list of names `self.selection` (if not given,
        # simply return all source functions).
        selection = group if not selection else selection
        for key in selection:
            S = group[key][()]

            container.sources[key] = S

        return container


class CLASSAdapter:
    def __init__(self, cosmo):
        self.cosmo = cosmo

    def get_inputs(self):
        container = InputContainer()

        cosmo = self.cosmo
        cosmo_params = cosmo.nn_cosmological_parameters()

        ###################################################################
        # 1) process thermodynamical + background quantities
        ###################################################################

        thermos = cosmo.get_thermos_for_NN()

        container.tau_th = thermos["tau"]

        thermos_params = cosmo.get_current_derived_parameters([
            "tau_rec",
            "z_reio",
            "rs_rec",
            "ds_rec",
            "ra_rec",
            "da_rec",
            "rd_rec"])
        tau_reio = thermos_params["tau_reio"] = cosmo.tau_of_z(thermos_params["z_reio"])
        del thermos_params["z_reio"]
        container.scalars.tau_rec = thermos_params["tau_rec"]

        container.scalars.update(thermos_params)

        container.scalars.tau_eisw_lisw_50 = cosmo.tau_of_z(50)
        container.scalars.tau_eisw_lisw_120 = cosmo.tau_of_z(120)

        z_bg = cosmo.get_bg_z()
        tau_bg = cosmo.get_bg_tau()

        tau_eq, a_eq, H_eq = cosmo.get_quantities_at_RM_equality()

        container.scalars.tau_eq = tau_eq
        container.scalars.a_eq = a_eq
        container.scalars.H_eq = H_eq

        container.scalars.rs_drag = self.cosmo.rs_drag_nn()

        bg_nn = cosmo.get_backgrounds_for_NN()

        tau_bg = bg_nn["tau"]
        r_s = bg_nn["r_s"]
        rho_b = bg_nn["rho_b"]
        rho_g = bg_nn["rho_g"]

        container.tau_bg = tau_bg
        container.background.r_s = r_s
        container.background.rho_b = rho_b
        container.background.rho_g = rho_g
        container.background.D = bg_nn["D"]

        container.thermo.e_kappa = thermos["e_kappa"]
        container.thermo.r_d = thermos["r_d"]

        container.thermo.g = thermos["g"]
        container.thermo.g_prime = thermos["g_prime"]
        container.thermo.g_reco = thermos["g_reco"]
        container.thermo.g_reco_prime = thermos["g_reco_prime"]
        container.thermo.g_reio = thermos["g_reio"]
        container.thermo.g_reio_prime = thermos["g_reio_prime"]

        container.scalars.k_min = self.cosmo.k_min()
        container.tau_source = self.cosmo.get_tau_source()

        ###################################################################
        # 2) process cosmological_parameters
        ###################################################################
        for key, value in cosmo_params.items():
            # Prefix cosmological parameters with 'cosmos/' to distinguish
            # the cosmological parameter "tau_reio" (i.e. the optical depth
            # to reionization) from the thermodynamical parameter "tau_reio"
            # (the conformal time of reionization)
            # This is done in order to be consistent with the naming scheme
            # of input parameters of CLASS
            container.cosmos["raw_cosmos/" + key] = value
            container.cosmos["cosmos/" + key] = value

        return container

    def get_outputs(self, selection=None):
        container = TargetContainer()
        sources, container.k_source, container.tau_source = self.cosmo.get_sources()

        selection = sources.keys() if not selection else selection
        for key in selection:
            container.sources[key] = sources[key]

        return container


