import numpy as np
import matplotlib.pyplot as plt 
import argparse 

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("history_file", type=str, help="File containing the training history (i.e. loss, val_loss, ...)")
    parser.add_argument("-l", "--loss-column", type=int, required=False, default=0, help="column number of column containining the (training) loss")
    parser.add_argument("-v", "--validation-loss-column", type=int, required=False, default=1, help="column number of column containining the validation loss")
    parser.add_argument("-g", "--show-grid", action="store_true", help="Plot a grid")
    parser.add_argument("--log", action="store_true", help="logarithmic y axis")

    args = parser.parse_args()
    data = np.genfromtxt(args.history_file)

    fig, ax = plt.subplots()

    ax.plot(data[:, args.loss_column], label="training loss", ls="-", marker="o")
    ax.plot(data[:, args.validation_loss_column], label="validation loss", ls="--", marker="o")

    if args.log:
        ax.set_yscale("log")

    ax.set_xlabel("Epoch")
    if args.show_grid:
        ax.grid()
    ax.legend()

    plt.show()

