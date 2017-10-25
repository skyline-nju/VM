"""
    Make a movie to show how the defect pair evolves with increasing time.
"""

import os
import numpy as np
import platform
import sys
import matplotlib
if platform.system() is "Windows":
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation
    from matplotlib.colors import hsv_to_rgb
    plt.rcParams['animation.ffmpeg_path'] = r"D:\ffmpeg\bin\ffmpeg"
else:
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation
    from matplotlib.colors import hsv_to_rgb
    dest = "/ffmpeg-3.3-64bit-static/ffmpeg"
    path1 = "/home-gk/users/nscc1185/Applications"
    path2 = "/home-yw/users/nsyw449_YK/dy/Program"
    if os.path.exists(path1):
        plt.rcParams['animation.ffmpeg_path'] = path1 + dest
    elif os.path.exists(path2):
        plt.rcParams['animation.ffmpeg_path'] = path2 + dest
    else:
        print("Error, cannot find ffmpeg")
        sys.exit()


def get_rgb(num, svx, svy, nrows, ncols, vmax=10):
    theta = np.arctan2(svy, svx)
    module = np.sqrt(svx**2 + svy**2)
    H = theta / (np.pi * 2) + 0.5
    max_module = module.max()
    if max_module < vmax:
        V = module / max_module
    else:
        V = module / vmax
        V[V > 1] = 1
    H = H.reshape(nrows, ncols)
    V = V.reshape(nrows, ncols)
    S = np.ones_like(H)
    HSV = np.dstack((H, S, V))
    RGB = hsv_to_rgb(HSV)
    return RGB


def set_dt(t):
    """ Set time interval between two frames. """
    if t < 30:
        dt = 1
    elif t < 60:
        dt = 2
    elif t < 120:
        dt = 4
    elif t < 240:
        dt = 8
    elif t < 480:
        dt = 16
    elif t < 960:
        dt = 32
    else:
        dt = 64
    return dt


def set_range(t, v0, nrows, ncols, l):
    """ Set the region to show. """
    d = 64 // l
    rowc = nrows // 2
    colc = ncols // 2
    if t < 100 or v0 == 0:
        row1, row2 = rowc - d, rowc + d
        col1, col2 = colc - d, colc + d
    elif t < 200:
        row1, row2 = rowc - d * 2, rowc + d * 2
        col1, col2 = colc - d * 2, colc + d * 2
    elif t < 800:
        row1, row2 = rowc - d * 4, rowc + d * 4
        col1, col2 = colc - d * 4, colc + d * 4
    elif t < 1600:
        row1, row2 = rowc - d * 8, rowc + d * 8
        col1, col2 = colc - d * 8, colc + d * 8
    else:
        row1, row2 = 0, nrows
        col1, col2 = 0, ncols
    box = [row1 * l, row2 * l, col1 * l, col2 * l]
    return col1, col2, row1, row2, box


def make_movie(L, rho0, d0, spin_num, nstep, eta, eps=0, seed=1, v0=0.5, l=1):
    """ Show the evolution of defect pairs.

    Left pannel for density fields while right pannel for orientational fields.

    Parameters:
    --------
        L : int
            System size.
        rho0 : float
            Particle number density.
        d0 : float
            Half of separation betwwen two defects.
        spin_num : int
            Rototion angle: 0, spin_num * PI / 2
        nstep : int
            Total time steps to run.
        eta : float
            Noise strength.
        eps : float, optional
            Strength of disorder.
        seed : int, optional
            Seed of random number generator.
        v0 : float, optional
            Velocity of particles.
        l : int, optional
            Boxes size for coarse-graining.
    """

    # create defect pair
    rot_angle = spin_num * np.pi / 2
    x, y, phi = defect.create_defect_pair(d0, rho0, L, rot_angle=rot_angle)
    num, svx, svy = defect.coarse_grain(x, y, phi, L, L, l)
    # initialize VMpy
    VMpy.ini(x, y, np.cos(phi), np.sin(phi), seed, v0, eta, eps, L, L)

    # initialize movie writter
    FFMpegWriter = animation.writers["ffmpeg"]
    writer = FFMpegWriter(fps=2, metadata=dict(artist="Matplotlib"))

    # set outfile name
    if spin_num == 0:
        mode = "a"
    elif spin_num == 1:
        mode = "b"
    elif spin_num == 2:
        mode == "c"
    else:
        mode = "d"
    filename = "data/L=%d_rho=%g_eta=%g_eps=%g_v=%g_%s.mp4" % (L, rho0, eta,
                                                               eps, v0, mode)
    fig = plt.figure(figsize=(12, 7))
    with writer.saving(fig, filename, dpi=100):
        # show inital configuration of defect pair
        for i in range(6):
            ax = plt.subplot(111)
            defect.show_defects(50, 30, d0, rot_angle, ax, L / 2, L / 2, True)
            title = "Inital state: system size $L=%d$, " % L \
                + "defect spacing $d_0=%d$, " % (2 * d0) \
                + "uniform density %d,\n" % (rho0) \
                + r"velocity $v_0=%g,\ $" % v0 + r"$\eta = %g,\ $" % eta \
                + r"$\epsilon=%g$" % eps
            ax.set_title(title, fontsize="xx-large")
            writer.grab_frame()
            plt.clf()

        # plot two panels
        ax1 = plt.subplot(121)
        ax2 = plt.subplot(122)
        im0 = ax1.imshow(
            np.zeros((L, L)),
            animated=True,
            extent=[0, L, 0, L],
            origin="lower",
            vmin=0,
            vmax=rho0 * 2.5)
        im1 = ax2.imshow(
            np.zeros((L, L)),
            animated=True,
            extent=[0, L, 0, L],
            origin="lower",
            vmin=-np.pi,
            vmax=np.pi,
            cmap="hsv")
        plt.colorbar(im0, ax=ax1, orientation="horizontal", extend="max")
        ticks = [-np.pi, -0.5 * np.pi, 0, 0.5 * np.pi, np.pi]
        labels = [r"$-\pi$", r"$-\pi/2$", r"$0$", r"$\pi/2$", r"$\pi$"]
        cb2 = plt.colorbar(im1, ax=ax2, orientation="horizontal", ticks=ticks)
        cb2.ax.set_xticklabels(labels)
        ax1.set_title("density")
        ax2.set_title("orientation")
        title = plt.suptitle("")
        title_template = r"$L=%d,\ \rho_0=%g, \eta=%g,\ \epsilon=%g,\ " \
            + r"v_0=%g,\ t=%d,\ \langle\phi\rangle=%.5f$"
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        ncols = L // l
        nrows = L // l
        t = 0
        while t <= nstep:
            print("t=%d" % t)
            dt = set_dt(t)
            if t > 0:
                VMpy.run(dt)
                VMpy.get_coarse_grained_snap(num, svx, svy, l)
            t += dt
            order_para = np.sqrt(np.sum(svx)**2 + np.sum(svy)**2) / (
                L**2 * rho0)
            col1, col2, row1, row2, box = set_range(t, v0, nrows, ncols, l)
            im0.set_data(num.reshape(nrows, ncols)[row1:row2, col1:col2])
            # phi2d = np.arctan2(svy, svx).reshape(nrows, ncols)
            # im1.set_data(phi2d[ymin:ymax, xmin:xmax])
            RGB = get_rgb(num, svx, svy, nrows, ncols, rho0 * 2.5)
            im1.set_data(RGB[row1:row2, col1:col2])
            title.set_text(title_template % (L, rho0, eta, eps, v0, t,
                                             order_para))
            im0.set_extent(box)
            im1.set_extent(box)
            writer.grab_frame()


if __name__ == "__main__":
    import defect
    import VMpy
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-N", "--nstep", type=int, help="Total time steps")
    parser.add_argument("--rho0", type=float, help="Particle density")
    parser.add_argument("-L", type=int, help="system size")
    parser.add_argument("-s", type=int, default=-1, help="Random number seed")
    parser.add_argument("-v", type=float, default=0.5, help="velocity")
    parser.add_argument(
        "--eta", type=float, default=0, help="strength of noise")
    parser.add_argument(
        "--eps", type=float, default=0, help="strength of disorder")
    parser.add_argument(
        "-l", type=int, default=1, help="boxes size for coarse-graining")
    parser.add_argument(
        "-a", "--angle", type=int, default=0, help="rotate angle: 0, 1, 2, 3")
    parser.add_argument(
        "-d", type=int, default=10, help="Half of defects spacing")
    arg = parser.parse_args()
    if arg.s < 0:
        import datetime
        seed = int(datetime.datetime.now().strftime("%m%d%H%M%S"))
    else:
        seed = arg.s
    make_movie(arg.L, arg.rho0, arg.d, arg.angle * np.pi / 2, arg.nstep,
               arg.eta, arg.eps, seed, arg.v, arg.l)
