'''
Calculate and show the order parameter for condensation
'''

import numpy as np
import matplotlib.pyplot as plt
from read_gsd import get_n_frames, read_frames

root_sohrab = "/run/user/1148/gvfs/sftp:host=sohrab003,user=yduan/scratch03.local/yduan"
root_tahmineh = "/run/user/1148/gvfs/sftp:host=tahmineh002,user=yduan/scratch03.local/yduan"

def cal_pos_order_para(fname):
    n_frames = get_n_frames(fname)
    t_arr = np.zeros(n_frames, int)
    psi_arr = np.zeros(n_frames)
    theta_arr = np.zeros(n_frames)

    frames = read_frames(fname)
    for i, frame in enumerate(frames):
        Lx = frame.configuration.box[0]
        t_arr[i] = frame.configuration.step
        x = frame.particles.position[:, 0]
        psi = x / Lx * np.pi * 2
        psi_arr[i] = np.sqrt(np.mean(np.cos(psi))**2 + np.mean(np.sin(psi))**2)
        theta = frame.particles.position[:, 2]
        theta_arr[i] = np.arctan2(np.mean(np.sin(theta)), np.mean(np.cos(theta)))
    return t_arr, psi_arr, theta_arr
    

def plot_psi_theta_varied_alpha(rho0):
    folder = f"{root_sohrab}/AsymVM/L512_128_half_self"

    Lx = 512
    Ly = 128
    eta = 0.03
    seed = 2001

    if rho0 == 3:
        a_arr = np.array([0, 0.01, 0.05, 0.1])
    elif rho0 == 0.65:
        a_arr = np.array([0.06])
    elif rho0 == 4:
        a_arr = np.array([0, 0.01, 0.05])
    elif rho0 == 5:
        a_arr = np.array([0, 0.01, 0.05])

    
    fig, (ax1, ax2) = plt.subplots(2, 1, constrained_layout=True, sharex=True)
    for a in a_arr:
        basename = f"L{Lx:d}_{Ly:d}_e{eta:g}_r{rho0:g}_a{a:g}_s{seed:d}.gsd"
        fname = f"{folder}/{basename}"
        t_arr, psi_arr, theta_arr = cal_pos_order_para(fname)

        if rho0 == 3 and a == 0.1:
            psi_arr = psi_arr[3:]
            theta_arr = theta_arr[3:]
            t_arr = t_arr[3:] - t_arr[3]
        if rho0 == 4 and a == 0.05:
            psi_arr = psi_arr[4:]
            theta_arr = theta_arr[4:]
            t_arr = t_arr[4:] - t_arr[4]
        if rho0 == 5 and a == 0.05:
            psi_arr = psi_arr[2:]
            theta_arr = theta_arr[2:]
            t_arr = t_arr[2:] - t_arr[2]
        ax1.plot(t_arr, psi_arr, "-o", label=r"$%g$" % a, ms=2)
        ax2.plot(t_arr, theta_arr, '-o', label=r"$%g$" % a, ms=2)

    ax1.legend(title=r"$\alpha=$")
    ax2.set_xlabel(r"$t$", fontsize="large")
    ax1.set_ylabel(r"$\psi$", fontsize="large")
    ax2.set_ylabel(r"$\theta_m$", fontsize="large")
    ax2.set_xlim(0, 1.8e6)
    ax1.set_ylim(ymin=0.9, ymax=1)
    ax2.axhline(0, linestyle="dotted", color="tab:grey")

    fig.suptitle(r"$\rho_0=%g,\eta=%g, L_x=%d, L_y=%d$" % (rho0, eta, Lx, Ly), fontsize="x-large")
    plt.show()
    plt.close()


if __name__ == "__main__":
    plot_psi_theta_varied_alpha(rho0=3)    
