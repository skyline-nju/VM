import numpy as np
import matplotlib.pyplot as plt
import os
import struct


def para2str(para: list, latex: bool=False) -> str:
    """ Transform the list of parameters into a formatted string. """
    eta, eps, rho, Lx, Ly, seed, t = para
    if latex:
        res = \
            r"$\eta=%g,\epsilon=%g,\rho_0=%g,L_x=%d,L_y=%d,\rm{seed}=%d,t=%d$"\
            % (eta, eps, rho, Lx, Ly, seed, t)
    else:
        res = "%g_%g_%g_%d_%d_%d_%08d" % (eta, eps, rho, Lx, Ly, seed, t)
    return res


def read(para: list) -> np.ndarray:
    """ Read x, y, theta from a binary file."""
    with open("s_%s.bin" % (para2str(para)), "rb") as f:
        buff = f.read()
        n = len(buff) // 12
        data = struct.unpack("%df" % (n * 3), buff)
        x, y, theta = np.array(data, dtype=np.float32).reshape(n, 3).T
    return x, y, theta


def write(para: list, coor: np.ndarray):
    """ Output data to a binary file. """
    data = coor.T
    file = "s_%s.bin" % (para2str(para))
    data.tofile(file)


def show_snap(para: list, **kwargs):
    """ Plot snapshot."""
    if "ax" in kwargs:
        ax = kwargs["ax"]
        is_show = False
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


def replicate_x(x0, y0, theta0, n, Lx, lim=None):
    """ Replicate snapshot along x direction.

        Parameters:
        --------
            x0: np.ndarray
                Original x coordinates
            y0: np.ndarray
                Original y coordiantes
            theta0: np.ndarray
                Original theta
            n: int
                Times to replicate
            Lx: int
                Original box size in x direction
            lim: list
                The left and right edges of band region

        Returns:
        --------
            x: np.ndarray
                New x coordinates
            y: np.ndarray
                New y coordinates
            theta: np.ndarray
                New theta
            Lx_new: int
                New box size in x direction
    """

    if lim is not None:
        # filter out the band region lim[0] < x < lim[1]
        mask = np.where((x0 < lim[0]) | (x0 > lim[1]))
        xc = x0[mask]
        yc = y0[mask]
        theta_c = theta0[mask]
        xc[xc < lim[0]] += (lim[1] - lim[0])
        dx = lim[0] + Lx - lim[1]
    else:
        xc = x0
        yc = y0
        theta_c = theta0
        dx = Lx

    N0 = x0.size
    Nc = xc.size
    N = N0 + n * Nc
    x, y, theta = np.zeros((3, N), dtype=np.float32)
    x[:N0] = x0
    y[:N0] = y0
    theta[:N0] = theta0
    Lx_new = Lx + n * dx
    for i in range(n):
        k1 = N0 + i * Nc
        k2 = N0 + (i + 1) * Nc
        x[k1:k2] = xc + (i + 1) * dx
        y[k1:k2] = yc
        theta[k1:k2] = theta_c
    return x, y, theta, Lx_new


def replicate(para, nx=0, ny=0, lim=None, reverse=False, coor=None):
    """ Replicat snapshot along x and y direction.

        Parameters:
        --------
            para: list
                List: eta, eps, rho, Lx, Ly, seed, t
            nx: int
                Times of replication along x direction
            ny: int
                Times of replication along y direction
            lim: list
                If not None, filter out the region between lim[0]
                and lim[1] when replicating.
            reverse: bool
                False for replicating along x direction first
            coor: np.ndarray
                Original coordination: x0, y0, theta0

        Returns:
        --------
            para_new: list
                eta, eps, rho_new, Lx_new, Ly_new, seed, t
            coor_new: np.ndarray
                Coordination of replicated snapshot

    """
    eta, eps, rho, Lx, Ly, seed, t = para
    if coor is not None:
        x0, y0, theta0 = coor
    else:
        x0, y0, theta0 = read(para)
    if not reverse:
        x1, y1, theta1, Lx_new = replicate_x(x0, y0, theta0, nx, Lx, lim)
        y2, x2, theta2, Ly_new = replicate_x(y1, x1, theta1, ny, Ly)
    else:
        y1, x1, theta1, Ly_new = replicate_x(y0, x0, theta0, ny, Ly, lim)
        x2, y2, theta2, Lx_new = replicate_x(x1, y1, theta1, nx, Lx)
    rho_new = x2.size / (Lx_new * Ly_new)
    para_new = [eta, eps, rho_new, Lx_new, Ly_new, seed, t]
    coor_new = np.array([x2, y2, theta2])
    write(para_new, coor_new)
    return para_new, coor_new


if __name__ == "__main__":
    os.chdir("VM/snap")
    eta = 0.35
    eps = 0
    rho = 1
    Lx = 150
    Ly = 100
    seed = 123
    t = 100000
    para = [eta, eps, rho, Lx, Ly, seed, t]
    show_snap(para)
