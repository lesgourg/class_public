import numpy as np
import torch

from classynet.data_providers import current_transformer
from classynet.training.dataset import SourceFileDataset
from classynet.models import ALL_NETWORK_CLASSES
from classynet.training.multi_trainer import MultiTrainer


class Trainer:
    def __init__(self, workspace):
        self.workspace = workspace

    def train_all_models(self, workers=0):
        return self.train_models(ALL_NETWORK_CLASSES, workers=workers)

    def train_model(self, model_class, workers=0):
        return self.train_models([model_class], workers=workers)

    def train_models(self, model_classes, workers=0):
        # TODO let user select device
        #device = torch.device("cuda")
        device = torch.device("cpu")

        k = self.workspace.loader().k()

        # convert k to torch tensor if not already
        if isinstance(k, np.ndarray):
            k = torch.from_numpy(k)

        k_net = k.to(device).float()

        # instantiate models
        models = [ctor(k_net) for ctor in model_classes]
        parameter_counts = {model.name(): sum(p.numel() for p in model.parameters()) for model in models}
        print("Network Parameter Counts:")
        print("\n".join("{}: {}".format(k, v) for k, v in parameter_counts.items()))
        print()

        trainer = MultiTrainer(self.workspace, models, device)

        in_transform, out_transform = current_transformer.get_pair(
            normalization_file=self.workspace.normalization_file,
            k=k
        )

        # Shared dataset arguments
        dataset_args = dict(
            input_transformer=in_transform,
            target_transformer=out_transform,
            input_selection=collect_model_inputs(models),
            source_function_selection=collect_model_outputs(models),
            k=k,
            # TODO optional limit?
        )

        training_dataset = SourceFileDataset(
            data_directory=self.workspace.training_data,
            **dataset_args
        )
        validation_dataset = SourceFileDataset(
            data_directory=self.workspace.validation_data,
            **dataset_args
        )
        test_dataset = SourceFileDataset(
            data_directory=self.workspace.test_data,
            **dataset_args
        )

        trainer.train(training_dataset, validation_dataset, test_dataset, loader_workers=workers)


def collect_model_inputs(models):
    """
    Returns the union of all required model inputs.
    """
    return set(inp for model in models for inp in model.required_inputs())

def collect_model_outputs(models):
    """
    Returns the union of all model targets.
    """
    return set(src for model in models for src in model.source_functions())
