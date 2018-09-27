"""
    Read and handle full serials of M.
    2018/9/23
"""
import numpy as np
import matplotlib.pyplot as plt
import glob
import os


def read_full_serials(filename):
    s = filename.replace(".dat", "").split("_")
    t_half = int(s[5])
    t_eq = int(s[6])

    with open(filename) as f:
        lines = f.readlines()[t_eq + 1:]
        m = np.array([float(line.replace("\n", "")) for line in lines])
        n = len(lines)
        n_period = n // (t_half * 2)
        Q = np.zeros(n_period)
        for i in range(n_period):
            Q[i] = np.mean(m[i * t_half * 2: (i + 1) * t_half * 2])
    return t_half, Q


def handle_full_serials(L, beta, eps, rho0, h0, seed):
    files = glob.glob("%d_%g_%g_%g_%g_*_2000_%d.dat" %
                      (L, beta, eps, rho0, h0, seed))
    t_half_arr = np.zeros(len(files))
    xi_arr = np.zeros_like(t_half_arr)
    for i, infile in enumerate(files):
        t_half_arr[i], Q = read_full_serials(infile)
        xi_arr[i] = L ** 2 * np.var(np.abs(Q[800:]))
        # plt.plot(np.arange(Q.size), Q)
        # plt.title(r"$t_{1/2}=%d$" % t_half_arr[i])
        # plt.show()
        # plt.close()
    plt.plot(t_half_arr, xi_arr, "o")
    plt.show()
    plt.close()


if __name__ == "__main__":
    os.chdir(r"D:\data\DPT\AIM\full_serials")
    L = 64
    beta = 1.9
    eps = 0.5
    rho0 = 6
    h0 = 0.1
    handle_full_serials(L, beta, eps, rho0, h0, 1)
