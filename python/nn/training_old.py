import os
from time import time
import json
from pprint import pprint

import numpy as np
import h5py as h5
import matplotlib.pyplot as plt

import torch

import k_standard
import tau_standard

import current_transformer
import dataset
import utils
import callbacks
import constants

class OutputConfig:
    def __init__(self, model_name, run_index, workspace_directory, create_plot=True):
        self.model_name = model_name
        self.run_index = run_index
        self.workspace_directory = workspace_directory

        self.create_plot = create_plot

        self.run_directory = os.path.join(self.workspace_directory, self.model_name, self.run_index)

        self.model_directory = os.path.join(self.run_directory, "models")
        self.model_checkpoint_path = os.path.join(self.model_directory, "model.checkpoint.{epoch}.h5")

        self.history_directory = os.path.join(self.run_directory, "history")
        self.history_log_directory = os.path.join(self.history_directory, "log")

        self.loss_file = os.path.join(self.history_directory, "loss.txt")
        self.learning_rate_file = os.path.join(self.history_directory, "lr.txt") 


class ModelTrainer:
    def __init__(self, model, training_config, training_generator, validation_generator, learning_rate_schedule=None, callbacks=None):
        self.model = model
        self.config = training_config
        self.training_generator = training_generator
        self.validation_generator = validation_generator
        self.learning_rate_schedule = learning_rate_schedule
        self.callbacks = [] if not callbacks else callbacks

    def train(self):
        config = self.config

        output_config = self.get_output_config()
        self.create_output_directories(output_config)

        history = self._train(output_config)

        self.save_model(self.model, output_config)
        if output_config.create_plot:
            self.plot_history(history, output_config)

        return history

    def get_output_config(self):
        config = self.config
        model_name = config["output"]["model_name"]
        run_index = config["output"]["run"]
        create_plot = config["output"].getboolean("plot", True)

        workspace_dir = os.path.expandvars("$CLASSNET_TRAINING_WORKSPACE")

        return OutputConfig(model_name, run_index, workspace_dir, create_plot)

    def create_output_directories(self, output_config):
        paths = [
                output_config.model_directory, 
                output_config.history_directory, 
                output_config.history_log_directory
                ]
        for path in paths:
            if not os.path.exists(path):
                print("Creating '{}' because it does not exist".format(path))
                os.makedirs(path)

    def save_model(self, model, output_config):
        final_model_path = os.path.join(output_config.model_directory, "model.h5")
        print("Saving model to {}.".format(final_model_path))
        model.save(final_model_path)

    def plot_history(self, history, output_config):
        plot_path = os.path.join(output_config.history_directory, "loss.png")
        import matplotlib as mpl
        mpl.use("Agg")
        plt.figure()
        plt.semilogy(history.history["loss"], label="loss", marker="o")
        plt.semilogy(history.history["val_loss"], label="validation loss", marker="x")
        plt.legend()
        plt.grid()
        plt.xlabel("Epoch")
        plt.ylabel("Loss")
        plt.savefig(plot_path, dpi=300, bbox_inches="tight")

    @property
    def epochs(self):
        return self.config["fit"].getint("epochs")

class TorchModelTrainer(ModelTrainer):
    def __init__(self, model, training_config):
        self.model = model
        self.config = training_config

    def _train(self, output_config):
        import net_fitter
        fitter = net_fitter.Fitter(self.model, use_cuda=True)

        training_dir = os.path.expandvars("$CLASSNET_TRAINING_DATA")
        validation_dir = os.path.expandvars("$CLASSNET_VALIDATION_DATA")

        training_limit = self.config["fit"].getint("training_limit", None)
        validation_limit = self.config["fit"].getint("validation_limit", None)

        dset = dataset.SourceFileDataset.for_model(self.model, constants.training_dir, 
                limit=training_limit)
        val_dset = dataset.SourceFileDataset.for_model(self.model, constants.validation_dir,
                limit=validation_limit)

        epochs = self.config["fit"].getint("epochs")
        num_workers = self.config["fit"].getint("workers", 0)

        cbacks = [
                callbacks.CSVLogger(output_config.loss_file),
                callbacks.CheckpointCallback(output_config.model_checkpoint_path)
                ]

        print(f"train for {epochs} epochs using {num_workers} workers.")
        print("training started...")

        time_start = time()
        history = fitter.fit(dset, val_dset, epochs=epochs, num_workers=num_workers, cbacks=cbacks)
        time_training = time() - time_start

        print("training finished.")
        print("training took {}s -> {}s / epoch.".format(
            time_training,
            time_training / epochs,
            ))

        # print(f"Writing training history to {output_config.loss_file}.")
        # history.as_dataframe().to_csv(output_config.loss_file, sep="\t")

        return history

    def save_model(self, model, output_config):
        final_model_path = os.path.join(output_config.model_directory, "model.pt")
        print("Saving model to {}.".format(final_model_path))
        torch.save(model.state_dict(), final_model_path)

    def _get_dataset(self, directory, input_transformer=None, target_transformer=None, limit=None):
        if input_transformer is None:
            input_transformer = current_transformer.get_input_transformer_normalizer()
        if target_transformer is None:
            target_transformer = current_transformer.get_target_transformer_normalizer(k=k_standard.K_STANDARD)

        return dataset.SourceFileDataset(
                directory, 
                input_transformer=input_transformer,
                target_transformer=target_transformer,
                source_function_selection=self.model.source_functions(), 
                input_selection=self.model.required_inputs(),
                k=k_standard.K_STANDARD,
                tau=self.model.tau_training(),
                limit=limit,
                )


class GeneratorModelTrainer(ModelTrainer):

    def get_worker_config(self):
        # Only use multiprocessing when files need to be read from disk 
        # (i.e. during full training).
        if not self.training_generator.is_cached:
            print("Not using caching; setting up multiprocessing for data feed.")
            if "workers" in self.config["fit"]:
                workers = self.config["fit"].getint("workers")
                print("Using {} workers (read from config)".format(workers))
            else:
                workers = int(os.environ.get("SLURM_NTASKS", os.cpu_count() // 2))

            return dict(
                use_multiprocessing=self.config["fit"].getboolean("use_multiprocessing", True),
                workers=workers,
            )
        else:
            return dict()

    def get_fit_arguments(self, output_config):
        config = self.config
        fit_kwargs = dict(
                generator=self.training_generator,
                validation_data=self.validation_generator,
                verbose=config["fit"].getint("verbose", 1),
                max_queue_size=config["fit"].getint("max_queue_size", 400),
                initial_epoch=config["fit"].getint("initial_epoch", 0),
                epochs=config["fit"].getint("epochs"),
        )

        fit_kwargs.update(self.get_worker_config())

        if not "callbacks" in fit_kwargs:
            fit_kwargs["callbacks"] = []

        fit_kwargs["callbacks"].extend([
                keras.callbacks.CSVLogger(output_config.loss_file, separator="\t"),
                LearningRateLogger(output_config.learning_rate_file), 
                # keras.callbacks.TensorBoard(history_log_dir),
                ])

        fit_kwargs["callbacks"].extend(self.callbacks)

        if config["output"].getboolean("checkpoints", True):
            fit_kwargs["callbacks"].append(
                keras.callbacks.ModelCheckpoint(output_config.model_checkpoint_path, verbose=0),
                )

        # LEARNING RATE SCHEDULER
        if self.learning_rate_schedule:
            print("Adding learning rate schedule: {}".format(self.learning_rate_schedule.__name__))
            lr_scheduler = keras.callbacks.LearningRateScheduler(self.learning_rate_schedule, verbose=1)
            fit_kwargs["callbacks"].append(lr_scheduler)

        return fit_kwargs

    def _train(self, output_config):
        fit_kwargs = self.get_fit_arguments(output_config)
        print("Training with following kwargs:")
        pprint(fit_kwargs)
        print("Starting training...")
        time_start = time()
        history = self.model.fit_generator(**fit_kwargs)
        time_training = time() - time_start
        print("finished training.")
        print("Model training took {}s -> {}s / epoch.".format(
            time_training,
            time_training / len(history.history["loss"])
            ))
        return history

