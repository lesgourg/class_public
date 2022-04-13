import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
from classynet.models import ALL_NETWORK_CLASSES

class HistoryPlotter:

    def __init__(self, workspace):
        self.workspace = workspace

    def plot_and_save(self):
        for cls in ALL_NETWORK_CLASSES:
            self.plot_training_history(cls.__name__)

    def plot_training_history(self, name):
        fig, ax = plt.subplots()

        epoch, loss, val_loss = np.genfromtxt(self.workspace.history_for(name), unpack=True)
        ln, = ax.semilogy(epoch + 1, loss, marker="o", label="Loss")
        ax.semilogy(epoch + 1, val_loss, marker="x", color=ln.get_color(), ls="--", label="Val. Loss")
        ax.set(xlabel="Epoch")
        ax.grid()
        ax.legend()

        fig_path = self.workspace.plots / "history_{}.png".format(name)
        print("Saving history to {}".format(fig_path))
        fig.savefig(fig_path, dpi=200, bbox_inches="tight")
        plt.close(fig)
        
        h5_path = self.workspace.plots / "history_{}.h5".format(name)
        print("Saving history as txt file to {}".format(h5_path))
        h5_file=h5.File(h5_path, "w")
        h5_file.create_dataset("epoch",data=epoch)
        h5_file.create_dataset("loss",data=loss)
        h5_file.create_dataset("val_loss",data=val_loss)
        h5_file.close()

