import os
from collections import namedtuple, defaultdict
import enum
import textwrap

import numpy as np
import torch
from tabulate import tabulate

from ..timer import Timer

class Phase(enum.Enum):
    Training = 0
    Validation = 1
    # Testing = 2

TrainingInfo = namedtuple("TrainingInfo", ["net", "optimizer", "criterion", "lr_scheduler"])
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

    def train(self, training_dataset, validation_dataset, loader_workers=0):
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

        # Find the net with the biggest number of epochs to train
        max_epochs = max(n.net.epochs() for n in self.nets)

        histories = [[] for _ in self.nets]

        for epoch in range(max_epochs):
            loss = self.run_epoch(Phase.Training, train_loader, epoch)
            val_loss = self.run_epoch(Phase.Validation, val_loader, epoch)

            # Log the losses to file
            for container, loss_col, val_loss_col, hist in zip(self.nets, loss.T, val_loss.T, histories):
                # Skip models that have finished training
                if epoch >= container.net.epochs():
                    continue
                # Save the updated history to file indicated by the
                # current workspace
                row = np.array([epoch, np.mean(loss_col), np.mean(val_loss_col)])
                hist.append(row)
                np.savetxt(self.workspace.history_for(container.net.name()), hist)

                # Also advance the LR scheduler
                container.lr_scheduler.step()

        # Finally, save the models.
        for cont in self.nets:
            net = cont.net
            path = self.workspace.model_path(net.name())
            print("Saving model {} to {}".format(net.name(), path))
            torch.save(net.state_dict(), path)


    def run_epoch(self, phase, loader, epoch, print_interval=20):
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
            with timer("multibatch"):
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
                for j, current in enumerate(self.nets):
                    # check if current model has finished training;
                    # if so, skip it
                    if epoch >= current.net.epochs():
                        continue
                    # Perform the forward pass through the network
                    with timer(("forward", current.net.name())):
                        y_net = current.net(x)
                        y_sliced = y_true[:, :, self.target_slices[j]]
                        losses[i, j] = loss = current.criterion(y_net, y_sliced)
                        running_loss[j] += loss.item()
                    # If we are training, perform backward pass + optimization step
                    if phase == Phase.Training:
                        with timer(("backward", current.net.name())):
                            current.optimizer.zero_grad()
                            loss.backward()
                            current.optimizer.step()

            if self.verbose and i % print_interval == print_interval - 1:
                running_loss /= print_interval
                self.print_info(phase, epoch, i, len(loader), running_loss, timer)
                running_loss = np.zeros(len(self.nets))

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
