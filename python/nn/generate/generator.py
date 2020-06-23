import json
from . import generate_data
from .generate_sources_files import generate_data, generate_parameters_and_data

class Generator:
    def __init__(self, workspace):
        self.workspace = workspace

    def generate_data(self, fixed, domain, training, validation, processes=None):
        self.write_manifest(fixed, domain.keys())

        generate_parameters_and_data(
            training,
            domain, fixed,
            self.workspace.training_data, processes=processes
        )
        generate_parameters_and_data(
            validation,
            domain, fixed,
            self.workspace.validation_data, processes=processes
        )

    def generate_data_for(self, fixed, training, validation, processes=None, fixed_training_only=None):
        """
        this method works with precomputed `domain` data.
        """

        self.write_manifest(fixed, training.keys())

        # this must happen AFTER writing the manifest!
        if fixed_training_only is not None:
            fixed = fixed.copy()
            fixed.update(fixed_training_only)

        generate_data(
            training,
            fixed,
            self.workspace.training_data, processes=processes
        )

        generate_data(
            validation,
            fixed,
            self.workspace.validation_data, processes=processes
        )

    def write_manifest(self, fixed, varying_names):
        # Save the manifest declaring the fixed and variable inputs
        # the data has been generated for
        with open(self.workspace.manifest, "w") as dest:
            json.dump(self.manifest(fixed, varying_names), dest)


    def manifest(self, fixed, varying_names):
        return {
            "fixed": fixed,
            "cosmological_parameters": list(varying_names)
        }
