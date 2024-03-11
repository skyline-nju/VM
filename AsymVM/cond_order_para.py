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
    c_arr = np.zeros(n_frames)

    frames = read_frames(fname)
    for i, frame in enumerate(frames):
        Lx = frame.configuration.box[0]
        t_arr[i] = frame.configuration.step
        x = frame.particles.position[:, 0]
        theta = x / Lx * np.pi * 2
        c_arr[i] = np.sqrt(np.mean(np.cos(theta))**2 + np.mean(np.sin(theta))**2)
    return t_arr, c_arr
    

if __name__ == "__main__":
    folder = f"{root_sohrab}/AsymVM/L512_128_half_self"

    Lx = 512
    Ly = 128
    eta = 0.03
    rho0 = 3
    seed = 2001
    a_arr = np.array([0, 0.01, 0.02, 0.05, 0.1])
    
    fig, ax = plt.subplots(1, 1, constrained_layout=True)
    for a in a_arr:
        basename = f"L{Lx:d}_{Ly:d}_e{eta:g}_r{rho0:g}_a{a:g}_s{seed:d}.gsd"
        fname = f"{folder}/{basename}"
        t_arr, c_arr = cal_pos_order_para(fname)

        if a == 0.1:
            c_arr = c_arr[3:]
            t_arr = t_arr[3:] - t_arr[3]
    
        ax.plot(t_arr, c_arr, "-o", label=r"$%g$" % a)

    ax.legend(title=r"$\alpha=$")
    ax.set_xlabel(r"$t$", fontsize="large")
    ax.set_ylabel(r"$\psi$", fontsize="large")
    fig.suptitle(r"$\rho_0=%g,\eta=%g, L_x=%d, L_y=%d$" % (rho0, eta, Lx, Ly), fontsize="x-large")
    plt.show()
    plt.close()
