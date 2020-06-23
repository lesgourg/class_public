from collections import defaultdict, OrderedDict
import time
import pandas as pd

import numpy as np

import torch

import callbacks

class History:
    def __init__(self):
        self.history = defaultdict(list)

    def as_dataframe(self):
        keys = ["loss", "val_loss"]
        epoch = np.arange(len(self.history["loss"]))
        df = pd.DataFrame.from_dict({"loss": self.history["loss"], "val_loss": self.history["val_loss"]})
        df.index.name = "epoch"
        return df

class Timer:
    def __init__(self):
        self.time_epoch = 0
        self.time_batch = 0
        self.time_loading = 0
        self.time_update = 0

        self.reset_running()

    def reset_running(self):
        self.time_batch_running = 0
        self.time_loading_running = 0
        self.time_update_running = 0

    def start_epoch(self):
        self._start_epoch = time.perf_counter()

    def end_epoch(self):
        self.time_epoch = time.perf_counter() - self._start_epoch

    def start_batch(self):
        self._start_batch = time.perf_counter()

    def end_batch(self):
        self.time_batch = time.perf_counter() - self._start_batch
        self.time_batch_running += self.time_batch

    def start_loading(self):
        self._start_loading = time.perf_counter()

    def end_loading(self):
        self.time_loading = time.perf_counter() - self._start_loading
        self.time_loading_running += self.time_loading

    def start_update(self):
        self._start_update = time.perf_counter()

    def end_update(self):
        self.time_update = time.perf_counter() - self._start_update
        self.time_update_running += self.time_update

class Fitter:

    def __init__(self, net, criterion=None, use_cuda=True):
        self.device = torch.device("cuda:0" if use_cuda else "cpu")
        self.net = net
        self.net.to(self.device)
        self.optimizer = self.net.optimizer()

        if not criterion:
            if hasattr(net, "criterion"):
                criterion = net.criterion()
                print("Using net.criterion():", criterion)
            else:
                criterion = torch.nn.MSELoss()

        self.output_order = self.net.source_functions()

        self.criterion = criterion

    def fit(self, dataset, val_dataset, epochs=1, print_interval=20, num_workers=0, cbacks=None):
        if cbacks is None:
            cbacks = []

        def collate(x):
            xs = [item[0] for item in x]
            ys = [item[1] for item in x]

            keys = xs[0].keys()

            if isinstance(ys[0], dict):
                ys = torch.cat([torch.stack([entry[key] for key in self.output_order], axis=2) for entry in ys])
            else:
                ys = torch.cat(ys, axis=0)

            xs = {k: torch.cat([item[k] for item in xs], axis=0) for k in keys}
            
            return xs, ys

        loader_kwargs = dict(num_workers=num_workers, batch_size=1, shuffle=True, pin_memory=True, collate_fn=collate)

        data_loader = torch.utils.data.DataLoader(dataset, **loader_kwargs)
        val_data_loader = torch.utils.data.DataLoader(val_dataset, **loader_kwargs)

        timer = Timer()
        history = History()

        if hasattr(self.net, "lr_scheduler"):
            lr_scheduler = self.net.lr_scheduler(self.optimizer)
            print("Model has LR schedule:", lr_scheduler)
        else:
            print("Not using LR schedule")
            lr_scheduler = None

        param_count = sum(p.numel() for p in self.net.parameters())
        print(f"Training model with {param_count} parameters.")

        state = callbacks.State(-1, self.net, self.optimizer, history)
        for cb in cbacks:
            cb.on_training_begin(state)

        for epoch in range(epochs):
            state = callbacks.State(epoch, self.net, self.optimizer, history)
            for cb in cbacks:
                cb.on_epoch_begin(state)

            timer.start_epoch()

            losses = self.training_phase(data_loader, epoch, print_interval)
            history.history["loss_full"].append(losses)
            history.history["loss"].append(np.mean(losses))

            timer.reset_running()

            val_losses = self.validation_phase(val_data_loader, epoch, print_interval)
            mean_val_loss = np.mean(val_losses)
            history.history["val_loss_full"].append(val_losses)
            history.history["val_loss"].append(mean_val_loss)

            timer.end_epoch()

            if lr_scheduler is not None:
                lr_scheduler.step()

            print("Epoch {}: initial loss = {:.3e}, final loss = {:.3e}, avg. loss = {:.3e}, avg. val. loss = {:.3e} ({:.2f}s)".format(
                epoch + 1,
                losses[0],
                losses[-1],
                np.mean(losses),
                mean_val_loss,
                timer.time_epoch
                ))

            state = callbacks.State(epoch, self.net, self.optimizer, history)
            for cb in cbacks:
                cb.on_epoch_end(state)

        return history

    def training_phase(self, data_loader, epoch, print_interval):
        timer = Timer()
        running_loss = 0
        timer.reset_running()
        data_loader_iter = iter(data_loader)

        losses = np.zeros(len(data_loader))
        self.net.train()
        for i in range(len(data_loader)):

            timer.start_batch()

            timer.start_loading()

            x, y_true = next(data_loader_iter)

            x = {k: v.to(self.device) for k, v in x.items()}
            y_true = y_true.to(self.device)

            timer.end_loading()

            timer.start_update()

            self.optimizer.zero_grad()
            y_net = self.net(x)
            losses[i] = loss = self.criterion(y_net, y_true)
            loss.backward()
            running_loss += loss.item()
            self.optimizer.step()

            timer.end_update()
            timer.end_batch()

            if i % print_interval == print_interval - 1:
                self.print_info("TRAINING", epoch, i, running_loss, timer, print_interval)
                timer.reset_running()
                running_loss = 0

        return losses
    
    def validation_phase(self, data_loader, epoch, print_interval):
        timer = Timer()
        val_losses = np.zeros(len(data_loader))
        running_loss = 0

        self.net.eval()

        with torch.no_grad():
            data_loader_iter = iter(data_loader)
            for i in range(len(data_loader)):

                timer.start_batch()

                timer.start_loading()

                x, y_true = next(data_loader_iter)

                x = {k: v.to(self.device) for k, v in x.items()}
                y_true = y_true.to(self.device)

                timer.end_loading()

                timer.start_update()

                y_net = self.net(x)
                val_losses[i] = loss = self.criterion(y_net, y_true)
                running_loss += loss.item()

                timer.end_update()

                timer.end_batch()

                if i % print_interval == print_interval - 1:
                    self.print_info("VALIDATION", epoch, i, running_loss, timer, print_interval)
                    timer.reset_running()
                    running_loss = 0

        return val_losses

    def print_info(self, phase, epoch, i, running_loss, timer, interval):
        print("{: <12}: [{}, {:5d}] loss: {:.3e}; {:.3f} ms ({:.3f} ms [loading] + {:.3f} ms [update])".format(
            phase,
            epoch + 1, i + 1, 
            running_loss / interval, 
            timer.time_batch_running * 1e3 / interval,
            timer.time_loading_running * 1e3 / interval,
            timer.time_update_running * 1e3 / interval,
            ))

