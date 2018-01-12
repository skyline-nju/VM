# calculate order parameters

import numpy as np
import matplotlib.pyplot as plt
import glob
import os


def read_phi(file, ncut=1500):
    print(file)
    with open(file) as f:
        lines = f.readlines()[ncut:]
        print(len(lines))
        phi = np.array([float(i.split("\t")[1]) for i in lines])
        return np.mean(phi)


if __name__ == "__main__":
    os.chdir(r"E:\data\random_torque\metric_free\phi")
    L = [32, 64, 128, 256]
    for l in L:
        files = glob.glob("%d_*.dat" % l)
        phi = np.zeros(len(files))
        eta = np.zeros_like(phi)
        for i, file in enumerate(files):
            eta[i] = float(file.split("_")[1])
            phi[i] = read_phi(file)
        phi = sorted(phi, key=lambda x: eta[np.where(phi == x)])
        eta = sorted(eta)
        plt.plot(eta, phi, "-o", label=r"$L=%d$" % l)
    plt.xlabel(r"$\eta$", fontsize="xx-large")
    plt.ylabel(r"$\langle \phi \rangle$", fontsize="xx-large")
    plt.legend(fontsize="xx-large")
    plt.tight_layout()
    plt.show()
    plt.close()
