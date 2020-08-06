import json
import numpy as np


def compose(*transformers):
    are_input_transformers = all(isinstance(t, InputTransformer) for t in transformers)
    are_target_transformers = all(isinstance(t, TargetTransformer) for t in transformers)
    if not (are_input_transformers ^ are_target_transformers):
        raise TypeError("transformers must _either_ be all InputTransformer instances _or_ all TargetTransformer instances.")

    if are_input_transformers:
        return CompositeInputTransformer(transformers)
    else:
        return CompositeTargetTransformer(transformers)


class Transformer:
    pass


class InputTransformer(Transformer):

    def transform_inputs(self, inputs, *args, **kwargs):
        return {key: self.transform_input(key, value, *args, **kwargs) for key, value in inputs.items()}

    def transform_input(self, name, input_, *args, **kwargs):
        return self.transform_inputs({name: input_}, *args, **kwargs)[name]


class TargetTransformer(Transformer):

    def transform_targets(self, targets, *args, **kwargs):
        return {key: self.transform_target(key, value, *args, **kwargs) for key, value in targets.items()}

    def untransform_targets(self, targets, *args, **kwargs):
        return {key: self.untransform_target(key, value, *args, **kwargs) for key, value in targets.items()}

    def transform_target(self, name, target, *args, **kwargs):
        return self.transform_targets({name: target}, *args, **kwargs)[name]

    def untransform_target(self, name, target, *args, **kwargs):
        return self.untransform_targets({name: target}, *args, **kwargs)[name]


class CanonicalInputTransformer(InputTransformer):

    def transform_inputs(self, inputs):
        inputs = inputs.copy()
        if "reference" in inputs:
            inputs["reference"] = self._transform_reference(inputs, inputs["reference"])
        return inputs

    def _transform_reference(self, inputs, reference):
        if inputs["k"].ndim == 2:
            k = inputs["k"][0]
        else:
            k = inputs["k"]

        reference = inputs["reference"][:, :, 0]
        reference_transformed = CanonicalTargetTransformer(k, inputs)._transform_phi_plus_psi(reference)
        return reference_transformed[:, :, None]


class RelativeInputTransformer(InputTransformer):

    TRANSFORM_QUANTITIES = (
            "cos", "sin", "j1", "j2", "reference",
            "g", "g_prime", "g_reco", "g_reco_prime", "g_reio", "g_reio_prime", "e_kappa"
            )

    def __init__(self, reference_inputs):
        self.reference_inputs = reference_inputs

    def transform_input(self, name, input_, *args, **kwargs):
        if self._should_transform(name):
            return input_ - self.reference_inputs[name]
        else:
            return input_

    def _should_transform(self, key):
        return key in RelativeInputTransformer.TRANSFORM_QUANTITIES

class RelativeTargetTransformer(TargetTransformer):

    def __init__(self, reference_targets):
        self.reference_targets = reference_targets

    def transform_targets(self, targets, *args, **kwargs):
        return {key: targets[key] - self.reference_targets[key] for key in targets}

    def untransform_targets(self, targets, *args, **kwargs):
        return {key: targets[key] + self.reference_targets[key] for key in targets}

class CompositeInputTransformer(InputTransformer):

    def __init__(self, transformers):
        assert all(isinstance(transformer, InputTransformer) for transformer in transformers)
        self.transformers = transformers

    def transform_inputs(self, inputs, *args, **kwargs):
        for transformer in self.transformers:
            inputs = transformer.transform_inputs(inputs, *args, **kwargs)
        return inputs


class CompositeTargetTransformer(TargetTransformer):

    def __init__(self, transformers):
        assert all(isinstance(transformer, TargetTransformer) for transformer in transformers)
        self.transformers = transformers

    def transform_targets(self, targets, *args, **kwargs):
        for transformer in self.transformers:
            targets = transformer.transform_targets(targets, *args, **kwargs)
        return targets

    def untransform_targets(self, targets, *args, **kwargs):
        for transformer in reversed(self.transformers):
            targets = transformer.untransform_targets(targets, *args, **kwargs)
        return targets


class CanonicalTargetTransformer(TargetTransformer):

    def __init__(self, k, inputs=None):
        self.k = k
        self.inputs = inputs

    def transform_targets(self, targets, inputs=None):
        if inputs is not None:
            self.inputs = inputs

        # TODO RE-ENABLE AFTER TESTING??
        # targets = targets.copy()
        # if "phi_plus_psi" in targets:
        #     pp = targets["phi_plus_psi"]
        #     targets["phi_plus_psi"] = self._transform_phi_plus_psi(pp)
        # if "phi" in targets:
        #     targets["phi"] = self._transform_phi_plus_psi(targets["phi"])
        # if "psi" in targets:
        #     targets["psi"] = self._transform_phi_plus_psi(targets["psi"])
        if "delta_m" in targets:
            delta_m = targets["delta_m"]
            targets["delta_m"] = self._transform_delta_m(delta_m)

        return targets

    def untransform_targets(self, targets, inputs=None):
        if inputs is not None:
            self.inputs = inputs

        # TODO RE-ENABLE AFTER TESTING??
        # targets = targets.copy()
        # if "phi_plus_psi" in targets:
        #     pp = targets["phi_plus_psi"]
        #     targets["phi_plus_psi"] = self._untransform_phi_plus_psi(pp)
        # if "phi" in targets:
        #     targets["phi"] = self._untransform_phi_plus_psi(targets["phi"])
        # if "psi" in targets:
        #     targets["psi"] = self._untransform_phi_plus_psi(targets["psi"])
        if "delta_m" in targets:
            delta_m = targets["delta_m"]
            targets["delta_m"] = self._untransform_delta_m(delta_m)

        return targets

    def _get_k_eq(self):
        a_eq = self.inputs["a_eq"][0]
        H_eq = self.inputs["H_eq"][0]
        return a_eq * H_eq

    def _get_k_th(self, k_eq):
        """
        :param k_eq: k of equality
        :return: k_threshold for normalization of delta_m and phi+psi
        """
        return 3.1 * k_eq

    def _transform_delta_m(self, value):
        k = self.k
        D = self.inputs["D"]
        k_eq = self._get_k_eq()

        return value / D[:, None]

        k_th = self._get_k_th(k_eq)
        mask = k < k_th

        scaled_factor = 1.8/(13.41*k_eq)
        result = value / D[:, np.newaxis]
        result[:, mask] *= (k_th / k[mask][np.newaxis, :])**2
        result[:, ~mask] *= np.log(np.e+scaled_factor*k_th) / np.log(np.e+scaled_factor*k[~mask][np.newaxis, :])
        result *= -1
        return result

    def _untransform_delta_m(self, value):
        k = self.k
        D = self.inputs["D"]

        return value * D[None, :]

        k_eq = self._get_k_eq()
        k_th = self._get_k_th(k_eq)
        mask = k < k_th

        scaled_factor = 1.8/(13.41*k_eq)
        result = -value
        result[~mask] /= np.log(np.e+scaled_factor*k_th) / np.log(np.e+scaled_factor*k[~mask][:, np.newaxis])
        result[mask, :] /= (k_th / k[mask][:, np.newaxis])**2
        result *= D[np.newaxis, :]
        return result

    def _transform_phi_plus_psi(self, value):
        k = self.k
        k_eq = self._get_k_eq()
        k_th = self._get_k_th(k_eq)

        mask = k < k_th
        scaled_factor = 1.8/(13.41*k_eq)
        result = np.copy(value)
        result[:, ~mask] *= np.log(np.e+scaled_factor*k_th) / np.log(np.e+scaled_factor*k[~mask][np.newaxis, :]) * (k[~mask][np.newaxis, :] / k_th)**2
        return result

    def _untransform_phi_plus_psi(self, value):
        k = self.k
        k_eq = self._get_k_eq()
        k_th = self._get_k_th(k_eq)

        mask = k < k_th

        scaled_factor = 1.8/(13.41*k_eq)
        result = np.copy(value)
        result[~mask] /= np.log(np.e+scaled_factor*k_th) / np.log(np.e+scaled_factor*k[~mask][:, np.newaxis]) * (k[~mask][:, np.newaxis] / k_th)**2
        return result


class Normalizer:

    def __init__(self, maxima, minima):
        self.maxima = maxima
        self.minima = minima

        self.abs_maxima = {k: max(abs(maxima[k]), abs(minima[k])) for k in maxima}

    @classmethod
    def from_path(cls, path, *args, **kwargs):
        with open(path) as infile:
            data = json.load(infile)
            return cls(data["max"], data["min"], *args, **kwargs)


class CanonicalInputNormalizer(Normalizer, InputTransformer):

    def __init__(self, maxima, minima):
        Normalizer.__init__(self, maxima, minima)

    def transform_inputs(self, inputs):
        return {key : self._normalize_input(key, value) for key, value in inputs.items()}

    def _normalize_input(self, key, value):
        if self._is_already_normalized(key) or self._should_pass_trough(key):
            return value
        elif self._should_normalize_by_maximum(key):
            return value / self.abs_maxima[key]
        elif self._is_cosmological_parameter(key):
            return self._normalize_cosmological_parameter(key, value)
        elif self._is_conformal_time(key):
            return np.log10(value)
        elif key.startswith("t0_reco_basis"):
            return self._normalize_t0_reco_basis(key, value)
        elif key in self.abs_maxima:
            return value / self.abs_maxima[key]
        else:
            raise ValueError("Do not know how to normalize input quantity {}!".format(key))

    def _normalize_t0_reco_basis(self, key, value):
        if key == "t0_reco_basis_with_g":
            assert "g_reco" in self.abs_maxima
            return value / self.abs_maxima["g_reco"]
        elif key == "t0_reco_basis_with_g_prime":
            assert "g_reco_prime" in self.abs_maxima
            return value / self.abs_maxima["g_reco_prime"]
        else:
            raise ValueError("Unknown: {}".format(key))

    def _should_pass_trough(self, key):
        return key in set((
            "k",
            "k_min",
            "r_s",
            "k_eq",
            "k_d",
            "rs_drag",
            "t0_reco_approx_6",
            "t0_reco_approx_7",
            "t0_reco_approx_8",
        )) or key.startswith("test_")

    def _is_already_normalized(self, key):
        return key.startswith("raw") or key in (
                "rho_b",
                "rho_g",
                "R",
                "D",
                "e_kappa",
                "damping",
                "cos",
                "sin",
                "cosj1",
                "sinj2",
                "reference",
                "D_prime",
                "j1_j2",
                "interp_param",
                "j1",
                "j2",
                "si_ci",
                "si",
                "ci",
                "t0_reco_basis",
                "phi_prime",
                "psi",
                "psi_minus_phi",
            "z",
            "a_eq",
            "k_eq",
            "H",
            "z_d",
                )

    def _should_normalize_by_maximum(self, key):
        return key in ("g", "g_reco", "g_reio", "g_prime", "g_reco_prime", "g_reio_prime",
                "t0_sw_approx_1", "t0_sw_approx_2",
                "t0_sw_approx_3", "t0_sw_approx_4",
                "t0_sw_approx_5", "t0_sw_approx_6",
                "t2_reco_approx_1", "t2_reco_approx_2",
                "t0_reco_approx_1", "t0_reio_approx_1",

                "t0_reco_approx_1",
                "t0_reco_approx_2",
                "t0_reco_approx_3",
                "t0_reco_approx_4",
                "t0_reco_approx_5",
                )

    def _is_conformal_time(self, key):
        return key in (
                "tau",
                "tau_rec",
                "tau_reio",
                "tau_eq",
                "tau_relative_to_reco",
                "tau_relative_to_reio",
                "tau_eisw_lisw_50",
                "tau_eisw_lisw_120",
                "tau_eisw_lisw_50_relative_to_rec",
                "tau_eisw_lisw_50_relative_to_reio",
                "tau_eisw_lisw_120_relative_to_rec",
                "tau_eisw_lisw_120_relative_to_reio",
                )

    def _is_cosmological_parameter(self, key):
        return key.startswith("cosmos/")

    def _normalize_cosmological_parameter(self, key, value):
        # _, key = key.split("/")
        diff = self.maxima[key] - self.minima[key]
        if diff < 1e-10:
            return value
        else:
            result = (value - self.minima[key]) / diff
        # shift from range (0, 1) to (-1, 1)
        result = 2 * result - 1
        return result

class AbsMaxNormalizer(Normalizer, TargetTransformer):

    def __init__(self, maxima, minima):
        Normalizer.__init__(self, maxima, minima)

    def transform_target(self, key, value, *args, **kwargs):
        return value / self.abs_maxima[key]

    def untransform_target(self, key, value, *args, **kwargs):
        # phi_prime has the same normalization as phi (because it has been fit
        # to match grad(phi, tau))
        if key == "phi_prime":
            key = "phi"
        return value * self.abs_maxima[key]
