import os
from collections import namedtuple, defaultdict
import enum
import textwrap
import copy

import numpy as np
import torch
from tabulate import tabulate

from classynet.tools.utils import Timer
from classynet.training import training_dashboard

class Phase(enum.Enum):
    Training = 0
    Validation = 1
    Test = 2

TrainingInfo = namedtuple("TrainingInfo", ["net", "optimizer", "criterion", "lr_scheduler"])

NET_NAME_DICT = {
    'Net_phi_plus_psi': 'psi_plus_phi',
    'Net_ST0_ISM': 't0_isw',
    'Net_ST0_Reco': 't0_reco_no_isw',
    'Net_ST0_Reio': 't0_reio_no_isw',
    'Net_ST1': 't1',
    'Net_ST2_Reco': 't2_reco',
    'Net_ST2_Reio': 't2_reio',
}

def build_info(net):
    optimizer = net.optimizer()
    lr_scheduler = net.lr_scheduler(optimizer)
    return TrainingInfo(net, optimizer, net.criterion(), lr_scheduler)


def get_collate(order):
    def collate(x):
        """
        Collation function for scenarios where batch_size != 1.
        """
        xs = [item[0] for item in x]
        ys = [item[1] for item in x]

        keys = xs[0].keys()

        if isinstance(ys[0], dict):
            ys = torch.cat([torch.stack([entry[key] for key in order], axis=2) for entry in ys])
        else:
            ys = torch.cat(ys, axis=0)

        xs = {k: torch.cat([item[k] for item in xs], axis=0) for k in keys}

        return xs, ys
    return collate

class MultiTrainer:

    def __init__(self, workspace, nets, device, verbose=True):
        self.workspace = workspace
        self.device = device
        self.verbose = verbose

        if not isinstance(nets, list) and not isinstance(nets, tuple):
            nets = [nets]

        # Make sure all nets are on the selected `device`
        for net in nets:
            net.to(self.device)

        self.nets = list(map(build_info, nets))

        # Since we combine the network training using a single DataLoader,
        # we will receive as training target a tensor of shape (#tau, #k, #source_functions),
        # so we need to keep track of:
        #   a) in which order the individual targets are put in this array
        #   b) which slices of this array belong to which network.
        # To do this, we keep track of the slices in the same order as the networks
        # in `self.nets` in `self.target_slices`, whose elements will be `slice`
        # instances or plain integer indices.
        # Also, we save the order in which the source functions are saved
        # into the tensor in `self.target_order`.
        self.target_slices = []
        self.target_order = []
        current_index = 0

        for container in self.nets:
            function_names = container.net.source_functions()
            self.target_order += function_names
            size = len(function_names)
            if size == 1:
                self.target_slices.append(current_index)
            else:
                self.target_slices.append(slice(current_index, current_index + size))
            current_index += len(function_names)

        # Find the net with the biggest number of epochs to train
        self.max_epochs = max(n.net.epochs() for n in self.nets)

        # Dict[str, List[float]]
        self.history = {c.net.name(): [] for c in self.nets}

    def train(self, training_dataset, validation_dataset, test_dataset, loader_workers=0):
        # arguments for the torch.utils.data.DataLoader instance
        print("Training with {} data loader workers.".format(loader_workers))
        loader_args = dict(
            num_workers=loader_workers,
            shuffle=True,
            # For fast GPU data transfer
            pin_memory=True,
            # The following two options are only needed if batch_size != 1
            batch_size=1,
            # TODO order?
            collate_fn=get_collate(self.target_order)
        )

        # Create the DataLoader instances from the given datasets
        train_loader = torch.utils.data.DataLoader(training_dataset, **loader_args)
        val_loader = torch.utils.data.DataLoader(validation_dataset, **loader_args)
        test_loader = torch.utils.data.DataLoader(test_dataset, **loader_args)

        histories = [[] for _ in self.nets]

        with training_dashboard.SimpleDashboard() as dashboard:
            self.dashboard = dashboard

            for epoch in range(self.max_epochs):
                loss = self.run_epoch(Phase.Training, train_loader, epoch)
                val_loss = self.run_epoch(Phase.Validation, val_loader, epoch)
                test_loss = self.run_epoch(Phase.Test, test_loader, epoch)

                # Log the losses to file
                for container, loss_col, val_loss_col, test_loss_col, hist in zip(self.nets, loss.T, val_loss.T, test_loss.T, histories):
                    name = container.net.name()
                    # Skip models that have finished training
                    if epoch >= container.net.epochs():
                        # Copy last validation loss
                        self.history[name].append(self.history[name][-1])
                        continue
                    self.history[name].append(np.mean(val_loss_col))
                    self.history[name].append(np.mean(test_loss_col))

                    # Save the updated history to file indicated by the
                    # current workspace
                    print(val_loss_col,test_loss_col)
                    row = np.array([epoch, np.mean(loss_col), np.mean(val_loss_col), np.mean(test_loss_col)])
                    hist.append(row)
                    np.savetxt(self.workspace.history_for(container.net.name()), hist)

                    # Also advance the LR scheduler
                    container.lr_scheduler.step()

                self.save_models(checkpoint=epoch)

        # Finally, save the models.
        self.save_models()

    def save_models(self, checkpoint=None):
        # here we store both the trained weights and an output normalization
        normalization_file = self.workspace.normalization_file()
        normalization = {key: np.max(abs(normalization_file['max'][key]),abs(normalization_file['min'][key])) for key in normalization_file['max'].keys()}

        for cont in self.nets:
            net = cont.net
            if checkpoint is None:
                path = self.workspace.model_path(net.name())
            else:
                path = self.workspace.model_path_checkpoint(net.name(), checkpoint)

            # we need to add the normalization constant to the networks
            my_state_dict = copy.deepcopy(net.state_dict())
            trans_net_name = NET_NAME_DICT[net.name()]
            my_state_dict['output_normalization'] = torch.tensor([normalization[trans_net_name]])
            print("Saving model {} to {}".format(net.name(), path))
            torch.save(my_state_dict, path)

    def run_epoch(self, phase, loader, epoch, print_interval=100):
        timer = Timer()

        # Prepare arrays keeping track of loss
        running_loss = np.zeros(len(self.nets))
        losses = np.zeros((len(loader), len(self.nets)))

        # Obtain an iterable over the data loader
        loader_iter = iter(loader)

        # Switch all nets into the respective mode, i.e. either training or evaluation
        for info in self.nets:
            if phase == Phase.Training:
                info.net.train()
            else:
                info.net.eval()

        for i in range(len(loader)):
            with timer("load data"):
                x, y_true = next(loader_iter)

                # If we are training a single network with a single input,
                # add a dummy channel axis so the slicing code below
                # doesn't break
                if y_true.ndim == 2:
                    y_true = y_true[:, :, None]

                x = {k: v.to(self.device) for k, v in x.items()}
                y_true = y_true.to(self.device)

            # iterate over networks
            stats = []
            for j, current in enumerate(self.nets):
                name = current.net.name()
                # check if current model has finished training;
                # if so, skip it
                if epoch >= current.net.epochs():
                    continue
                stat_dict = {"network": name}
                # Perform the forward pass through the network
                with timer(("forward", name)):
                    y_net = current.net(x)
                    y_sliced = y_true[:, :, self.target_slices[j]]
                    losses[i, j] = loss = current.criterion(y_net, y_sliced)
                    stat_dict["loss"] = loss.item()
                    running_loss[j] += loss.item()
                stat_dict["forward"] = timer.times[("forward", name)]
                # If we are training, perform backward pass + optimization step
                if phase == Phase.Training:
                    with timer(("backward", name)):
                        current.optimizer.zero_grad()
                        loss.backward()
                        current.optimizer.step()
                    stat_dict["backward"] = timer.times[("backward", name)]
                else:
                    stat_dict["backward"] = 0
                stats.append(stat_dict)

            if self.dashboard and i % print_interval == 0:
                self.dashboard.update(
                    phase,
                    epoch, self.max_epochs,
                    i, len(loader),
                    loading_time=timer.times["load data"],
                    live=stats,
                    history=self.history
                )

            # if self.verbose and i % print_interval == print_interval - 1:
            #     running_loss /= print_interval
            #     self.print_info(phase, epoch, i, len(loader), running_loss, timer)
            #     running_loss = np.zeros(len(self.nets))

        return losses

    def print_info(self, phase, epoch, i, batches, running_loss, timer):
        epoch_info = f"Epoch {epoch + 1: 3d}: {phase: <12}"
        batch_info = f"\tBatch {i} / {batches}"
        loading_info = "\tData Loading: {:.1f}ms".format(timer.times["load data"] * 1e3)

        rows = []
        for i, info in enumerate(self.nets):
            name = info.net.name()
            if epoch >= info.net.epochs():
                row = [name, "FINISHED", "0", "0"]
            else:
                f_time = timer.times[("forward", name)] * 1e3
                b_time = timer.times.get(("backward", name), 0) * 1e3
                row = [
                    name,
                    "{:.3e}".format(running_loss[i]),
                    round(f_time, 2),
                    round(b_time, 2),
                ]
            rows.append(row)

        table = tabulate(
            rows,
            headers=["model", "loss", "forward time (ms)", "backward time (ms)"],
            tablefmt="psql"
        )

        # os.system("clear")
        print(epoch_info)
        print("=" * len(epoch_info))
        print(batch_info)
        print(loading_info)
        print(textwrap.indent(table, "\t"))
        print()
