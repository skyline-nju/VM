"""
    Read and handle serials of period-averaged order parameters, Q.
    2018/9/23
"""
import numpy as np
import matplotlib.pyplot as plt
import glob
import os


def read_serials(filename, ncut=1000):
    s = filename.replace(".dat", "").split("_")
    t_half = int(s[5])

    with open(filename) as f:
        lines = f.readlines()[ncut:]
        Q = np.array([float(i.replace("\n", "")) for i in lines])
    return t_half, Q


def handle_all(L, beta, eps, rho0, h0, seed, axes=None):
    if axes is None:
        fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, nrows=1, figsize=(12, 4))
        flag_show = True
    else:
        ax1, ax2, ax3 = axes
        flag_show = False
    files = glob.glob("%d_%g_%g_%g_%g_*_%d.dat" %
                      (L, beta, eps, rho0, h0, seed))
    periods = np.zeros(len(files))
    qm = np.zeros_like(periods)
    xi = np.zeros_like(periods)
    Ul = np.zeros_like(periods)
    for i, file_i in enumerate(files):
        t_half, Q = read_serials(file_i)
        periods[i] = t_half * 2
        qm[i] = np.abs(np.mean(Q))
        xi[i] = L ** 2 * np.var(np.abs(Q))
        Ul[i] = 1 - np.mean(Q**4) / np.mean(Q**2) ** 2 / 3
    sort_idx = periods.argsort()
    periods = periods[sort_idx]
    qm = qm[sort_idx]
    xi = xi[sort_idx]
    Ul = Ul[sort_idx]
    ax1.plot(periods / 2, qm, "-o", label=r"$L=%d$" % L)
    ax2.plot(periods / 2, xi, "-o")
    ax3.plot(periods / 2, Ul, "-o")
    ax2.set_yscale("log")
    ax1.set_xlabel(r"$t_{1/2}$", fontsize="large")
    ax2.set_xlabel(r"$t_{1/2}$", fontsize="large")
    ax3.set_xlabel(r"$t_{1/2}$", fontsize="large")
    ax1.set_ylabel(r"$Q$", fontsize="large")
    ax2.set_ylabel(r"$\chi$", fontsize="large")
    ax3.set_ylabel(r"$U_L$", fontsize="large")
    if flag_show:
        plt.tight_layout()
        plt.show()
        plt.close()


if __name__ == "__main__":
    os.chdir(r"D:\data\DPT\AIM\q_serials")
    L = [64, 90, 128]
    beta = 1.9
    eps = 0.5
    rho0 = 6
    h0 = 0.1
    fig, axes = plt.subplots(ncols=3, nrows=1, figsize=(12, 4))
    for l in L:
        handle_all(l, beta, eps, rho0, h0, 1, axes)
    axes[0].legend()
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.suptitle(r"$\beta=%g,\ \epsilon=%g,\ \rho_0=%g, h_0=0.1$" %
                 (beta, eps, rho0), fontsize="xx-large")

    plt.show()
    plt.close()
