import os
import numpy as np
import matplotlib.pyplot as plt
from classynet.models import ALL_NETWORK_CLASSES
import os

#from class.python.nn.workspace import Workspace

class HistoryPlotter:

    def __init__(self, workspace):
        self.workspace = workspace

    def plot_and_save(self):
        for cls in ALL_NETWORK_CLASSES:
            if hasattr(self.workspace,'generations'):
                if cls.__name__ in self.workspace.generations.keys():
                    self.plot_training_history(cls.__name__)
            else:
                self.plot_training_history(cls.__name__)

    def plot_training_history(self, name):
        fig, ax = plt.subplots()
            
        epoch, loss, val_loss, test_loss = np.genfromtxt(self.workspace.history_for(name), unpack=True)
        ln, = ax.semilogy(epoch + 1, loss, marker="x", label="Loss", color='blue', ls="--")
        ax.semilogy(epoch + 1, val_loss, marker="x", color='green', ls="--", label="Val. Loss")
        ax.semilogy(epoch + 1, test_loss, marker="x", color='red', ls="--", label="Test Loss")
        ax.set(xlabel="Epoch")
        ax.grid()
        ax.legend()

        os.makedirs(self.workspace.plots / "training_history",exist_ok=True)
        if hasattr(self.workspace,'generations'):
            fig_path = self.workspace.plots / "training_history" / "history_{}_{}.png".format(name,self.workspace.generations[name])
        else:
            fig_path = self.workspace.plots / "training_history" / "history_{}.png".format(name)
        print("Saving history to {}".format(fig_path))
        fig.savefig(fig_path, dpi=200, bbox_inches="tight")
        plt.close(fig)
        
