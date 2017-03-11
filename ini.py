import numpy as np
import matplotlib.pyplot as plt
import os
import struct


def para2str(para, latex=False):
    eta, eps, rho, Lx, Ly, seed, t = para
    if latex:
        res = \
            r"$\eta=%g,\epsilon=%g,\rho_0=%g,L_x=%d,L_y=%d,\rm{seed}=%d,t=%d$"\
            % (eta, eps, rho, Lx, Ly, seed, t)
    else:
        res = "%g_%g_%g_%d_%d_%d_%08d" % (eta, eps, rho, Lx, Ly, seed, t)
    return res


def read(para):
    with open("s_%s.bin" % (para2str(para)), "rb") as f:
        buff = f.read()
        n = len(buff) // 12
        data = struct.unpack("%df" % (n * 3), buff)
        x, y, theta = np.array(data, dtype=np.float32).reshape(n, 3).T
    return x, y, theta


def write(para, coor):
    data = coor.T
    file = "s_%s.bin" % (para2str(para))
    data.tofile(file)


def show_snap(para, **kwargs):
    is_show = False
    if "ax" in kwargs:
        ax = kwargs["ax"]
    else:
        ax = plt.subplot(111)
        is_show = True
    if "coor" in kwargs:
        x, y, theta = kwargs["coor"]
    else:
        x, y, theta = read(para)
    ax.scatter(x, y, s=1, c=theta, cmap="hsv")
    ax.set_xlim(0, para[3])
    ax.set_ylim(0, para[4])
    ax.set_title(para2str(para, latex=True))
    if is_show:
        plt.show()
        plt.close()


def replicate(para, **kwargs):
    eta, eps, rho, Lx, Ly, seed, t = para
    if "coor" in kwargs:
        x0, y0, theta0 = kwargs["coor"]
    else:
        x0, y0, theta0 = read(para)

    full_copy = False
    if "xlim" in kwargs:
        xlim = kwargs["xlim"]
        ylim = [0, Ly]
        nx = kwargs["nx"]
        ny = 1
    elif "ylim" in kwargs:
        ylim = kwargs["ylim"]
        xlim = [0, Lx]
        ny = kwargs["ny"]
        nx = 1
    else:
        full_copy = True
        xlim = [0, Lx]
        ylim = [0, Ly]
        if "nx" in kwargs:
            nx = kwargs["nx"]
        else:
            nx = 1
        if "ny" in kwargs:
            ny = kwargs["ny"]
        else:
            ny = 1

    N0 = x0.size
    if full_copy:
        N_new = x0.size * nx * ny
        x_new = np.zeros(N_new, np.float32)
        y_new = np.zeros(N_new, np.float32)
        theta_new = np.zeros(N_new, np.float32)
        for j in range(ny):
            for i in range(nx):
                left = (i + j * nx) * N0
                right = (i + 1 + j * nx) * N0
                x_new[left:right] = x0 + i * Lx
                y_new[left:right] = y0 + j * Ly
                theta_new[left:right] = theta0
        para_new = [eta, eps, rho, Lx * nx, Ly * ny, seed, t]
        coor_new = np.array([x_new, y_new, theta_new])
    # elif ny == 1:
    #     mask = np.where((x0 < xlim[0]) & (x0 > xlim[1]))
    #     xc, yc, thetac = x0[mask], y0[mask], theta0[mask]
    #     xc[xc < xlim[0]] += (xlim[1] - xlim[0])
    #     lx = xlim[0] + Lx - xlim[1]
    #     N_new = N0 + (nx - 1) * xc.size
    #     x_new = np.zeros(N_new, np.float32)
    #     y_new = np.zeros(N_new, np.float32)
    #     theta_new = np.zeros(N_new, np.float32)
    #     for i in range(nx):
    #         left = (i)
    write(para_new, coor_new)
    return para_new


if __name__ == "__main__":
    os.chdir("snap")
    eta = 0.38
    eps = 0
    rho = 1
    Lx = 200
    Ly = 200
    seed = 123
    t = 100000
    para = [eta, eps, rho, Lx, Ly, seed, t]
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2)
    x, y, theta = read(para)
    coor = np.array([x, y, theta])
    nx = 2
    ny = 1
    show_snap(para, coor=coor, ax=ax1)
    para1 = replicate(para, nx=nx, ny=ny)
    show_snap(para1, ax=ax2)
    plt.show()
    plt.close()
