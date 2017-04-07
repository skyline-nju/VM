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


def read_coarse_grain(eta, eps, Lx, Ly, ncols, nrows, N, dt, seed, ff):
    """ Read coarse-grained snapshot. """

    file = "coarse/c%s_%g_%g_%d_%d_%d_%d_%d_%d_%d.bin" % (
        ff, eta, eps, Lx, Ly, ncols, nrows, N, dt, seed)
    ncell = ncols * nrows
    if ff == "Bbb":
        BLOCK_SIZE = 20 + ncell * 3
        format_str = "i2d%dB%db" % (N, 2 * N)
    elif ff == "iff":
        BLOCK_SIZE = 20 + ncell * 3 * 4
        format_str = "i2d%di%df" % (N, 2 * N)
    with open(file, "rb") as f:
        while True:
            block = f.read(BLOCK_SIZE)
            if not block:
                break
            else:
                data = struct.unpack(format_str, block)
                t = data[0]
                vxm = data[1]
                vym = data[2]
                num_ij = np.array(data[3:3 + ncell])
                vx_ij = np.array(data[3 + ncell:3 + 2 * ncell])
                vy_ij = np.array(data[3 + 2 * ncell:3 + 3 * ncell])
                if ff == "Bbb":
                    vx_ij /= 128
                    vy_ij /= 128
                frame = np.array([
                    t, vxm, vym, num_ij.reshape(nrows, ncols), vx_ij.reshape(
                        nrows, ncols), vy_ij.reshape(nrows, ncols)
                ])
                yield frame


def read_frame(eta, eps, Lx, Ly, N, dt, seed, ncols=None, nrows=None, ff=None):
    """ Read a frame from a large binary file."""
    if ncols is None:
        file = "snap_one/so_%g_%g_%d_%d_%d_%d_%d.bin" % (eta, eps, Lx, Ly, N,
                                                         dt, seed)
        BLOCK_SIZE = N * 3 * 4
        with open(file, "rb") as f:
            while True:
                block = f.read(BLOCK_SIZE)
                if not block:
                    break
                else:
                    frame = struct.unpack("%df" % (N * 3), block)
                    frame = np.array(frame, dtype=np.float32).reshape(N, 3).T
                    yield frame
    elif ff == "iff":
        file = "coarse/ciff_%g_%g_%d_%d_%d_%d_%d_%d_%d.bin" % (
            eta, eps, Lx, Ly, ncols, nrows, N, dt, seed)
        BLOCK_SIZE = N * 4
        with open(file, "rb") as f:
            while True:
                block = f.read(BLOCK_SIZE)
                if not block:
                    break
                else:
                    sum_num = np.array(
                        struct.unpack("%di" % N, block), dtype=int)
                block = f.read(BLOCK_SIZE * 2)
                if not block:
                    break
                else:
                    sum_vx, sum_vy = np.array(
                        struct.unpack("%df" % (2 * N), block),
                        dtype=np.float32).reshape(2, N)
                frame = np.array([
                    sum_num.reshape(nrows, ncols),
                    sum_vx.reshape(nrows, ncols), sum_vy.reshape(nrows, ncols)
                ])
                yield frame
    elif ff == "Bbb":
        file = "coarse/cBbb_%g_%g_%d_%d_%d_%d_%d_%d_%d.bin" % (
            eta, eps, Lx, Ly, ncols, nrows, N, dt, seed)
        BLOCK_SIZE = N
        with open(file, "rb") as f:
            while True:
                block = f.read(BLOCK_SIZE)
                if not block:
                    break
                else:
                    sum_num = np.array(
                        struct.unpack("%dB" % N, block), dtype=int)
                block = f.read(BLOCK_SIZE * 2)
                if not block:
                    break
                else:
                    sum_vx, sum_vy = np.array(
                        struct.unpack("%db" % (2 * N), block),
                        dtype=np.float32).reshape(2, N)
                sum_vx /= 128
                sum_vy /= 128
                frame = np.array([
                    sum_num.reshape(nrows, ncols),
                    sum_vx.reshape(nrows, ncols), sum_vy.reshape(nrows, ncols)
                ])
                yield frame


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
    # os.chdir("VM/snap")
    # eta = 0.35
    # eps = 0
    # rho = 1
    # Lx = 140
    # Ly = 200
    # seed = 33
    # t = 100000
    # para = [eta, eps, rho, Lx, Ly, seed, t]
    # show_snap(para)

    # os.chdir("VM")
    # eta = 0.35
    # eps = 0
    # rho = 1
    # Lx = 140
    # Ly = 200
    # N = 28000
    # dt = 10000
    # seed = 33
    # f1 = read_frame(
    #     eta, eps, Lx, Ly, N, dt, seed, ncols=Lx, nrows=Ly, ff="iff")
    # f2 = read_frame(
    #     eta, eps, Lx, Ly, N, dt, seed, ncols=Lx, nrows=Ly, ff="Bbb")

    # for frame1, frame2 in zip(f1, f2):
    #     n1, vx1, vy1 = frame1
    #     n2, vx2, vy2 = frame2
    #     print("vx error = %f" % (np.std(vx1 - vx2)))
    #     print("vy error = %f" % (np.std(vy1 - vy2)))
    #     plt.subplot(121)
    #     plt.contourf(vy1)
    #     plt.colorbar()
    #     plt.subplot(122)
    #     plt.contourf(vy2)
    #     plt.colorbar()
    #     plt.show()
    #     plt.close()


