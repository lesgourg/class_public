import sys
import numpy as np
import matplotlib.pyplot as plt


def main(argv):
    file = 'output/test_cl.dat'
    if len(argv) > 0:
        file = argv[0]
    dat = np.loadtxt(file).T
    ell = dat[0]
    c_ell = dat[1:]
 
    plt.figure(file)
    plt.xlabel(r'$\ell$')
    plt.ylabel(r'$\ell (\ell + 1) C_\ell$')
    plt.xscale('log')
    plt.yscale('log')
    for i,y in enumerate(c_ell):
        plt.plot(ell, y, label='%i'%i)
    plt.legend()
    plt.show()


if __name__ == "__main__":
    main(sys.argv[1:])