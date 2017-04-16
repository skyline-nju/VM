""" Prepare the snapshot for simulation.

    Generate new snapshot by coping and pasting a rectangle region of the
    original snapshot along x or y direction. Then initialize the similation
    from the new snapshot.
"""
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


def get_time_step(t_end, exponent, show=False):
    ts = [1]
    t_cur = 1
    while True:
        t_cur = t_cur * exponent
        if t_cur > t_end:
            break
        t_cur_int = int(np.round(t_cur))
        if t_cur_int > ts[-1]:
            ts.append(t_cur_int)
    if show:
        plt.subplot(121)
        plt.plot(ts, "-o")
        plt.subplot(122)
        plt.plot(ts, "-o")
        plt.yscale("log")
        plt.suptitle("n = %d" % (len(ts)))
        plt.show()
        plt.close()
    return ts


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


def replicate_x(x0, y0, theta0, Lx0, n, lim=None):
    """ Replicate snapshot along x direction.

        Parameters:
        --------
        x0, y0, theta0: np.ndarray
            Coordination and velocity angle of particles.
        Lx0: double
            System size in the x direction.
        n: int
            Times of replicating along x direction.
        lim: list
            The region to copy.

        Returns:
        --------
        x, y, theta: np.ndarray
            Coordiantion and velocity angle of new snapshot.
        Lx: double
            New system size in the x direction.
    """
    if lim is not None:
        xmin, xmax = lim
    else:
        xmin, xmax = 0, Lx0
    if xmin < xmax:
        x0 += (Lx0 - xmax)
        xmin += (Lx0 - xmax)
        x0[x0 >= Lx0] -= Lx0
    else:
        x0 -= xmax
        xmin -= xmax
        x0[x0 < 0] += Lx0
    mask = x0 >= xmin
    xc = x0[mask]
    yc = y0[mask]
    theta_c = theta0[mask]
    Lc = Lx0 - xmin
    x, y, theta = np.zeros((3, x0.size + n * xc.size), dtype=np.float32)
    x[:x0.size] = x0
    y[:x0.size] = y0
    theta[:x0.size] = theta0
    Lx = Lx0 + n * Lc
    for i in range(n):
        k1 = x0.size + i * xc.size
        k2 = x0.size + (i + 1) * xc.size
        x[k1:k2] = xc + (i + 1) * Lc
        y[k1:k2] = yc
        theta[k1:k2] = theta_c
    return x, y, theta, Lx


def replicate(para, nx=0, ny=0, lim=None, reverse=False, coor=None, out=False):
    """ Replicate snapshot along x and y direction.

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
    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1)
    show_snap(para, coor=[x0, y0, theta0], ax=ax1)
    if not reverse:
        x1, y1, theta1, Lx_new = replicate_x(x0, y0, theta0, Lx, nx, lim)
        y2, x2, theta2, Ly_new = replicate_x(y1, x1, theta1, Ly, ny)
    else:
        y1, x1, theta1, Ly_new = replicate_x(y0, x0, theta0, Ly, ny, lim)
        x2, y2, theta2, Lx_new = replicate_x(x1, y1, theta1, Lx, nx)
    rho_new = x2.size / (Lx_new * Ly_new)
    para_new = [eta, eps, rho_new, Lx_new, Ly_new, seed, t]
    coor_new = np.array([x2, y2, theta2])
    if out:
        write(para_new, coor_new)
    show_snap(para_new, coor=coor_new, ax=ax2)
    plt.show()
    plt.close()


if __name__ == "__main__":
    os.chdir("D:\\tmp")

    t = 8 * 100000
    para_list = [0.35, 0.02, 1, 220, 100, 123, t]
    replicate(para_list, 2, 0, lim=[110, 210], out=False)
    # get_time_step(500000, 1.06, show=True)
