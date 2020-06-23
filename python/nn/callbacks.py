import torch

class State:
    def __init__(self, epoch, model, optimizer, history):
        self.epoch = epoch
        self.model = model
        self.optimizer = optimizer
        self.history = history

class Callback:

    def on_training_begin(self, state):
        pass

    def on_epoch_begin(self, state):
        pass

    def on_epoch_end(self, state):
        pass

class CSVLogger(Callback):

    def __init__(self, output_path):
        super().__init__()

        self.output_path = output_path

    def on_epoch_end(self, state):
        print(f"CSVLogger: Epoch {state.epoch}, writing losses to {self.output_path}.")
        df = state.history.as_dataframe()
        df.to_csv(self.output_path, sep="\t")

class CheckpointCallback(Callback):
    def __init__(self, destination="model.checkpoint.{epoch}.pt"):
        self.destination = destination

    def on_epoch_end(self, state):
        path = self.destination.format(epoch=state.epoch)
        torch.save(state.model.state_dict(), path)

