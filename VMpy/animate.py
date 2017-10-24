"""
    Make a movie to show how the defect pair evolves with increasing time.
"""

import os
import numpy as np
import platform
import sys
import matplotlib
matplotlib.use("Agg")
if platform.system() is "Windows":
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation
    plt.rcParams['animation.ffmpeg_path'] = r"D:\ffmpeg\bin\ffmpeg"
else:
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation
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


def set_dt(t):
    """ Set time interval between two frames. """
    if t < 20:
        dt = 1
    elif t < 40:
        dt = 2
    elif t < 80:
        dt = 4
    elif t < 160:
        dt = 8
    elif t < 320:
        dt = 16
    elif t < 640:
        dt = 32
    else:
        dt = 64
    return dt


def make_movie(L, rho0, d0, rot_angle, nstep, eta, eps=0, seed=1, v0=0.5, l=1):
    """ Show the evolution of defect pairs.

    Parameters:
    --------
        L : int
            System size.
        rho0 : float
            Particle number density.
        d0 : float
            Half of separation betwwen two defects.
        rot_angle : float
            Rototion angle: 0, 1/2 PI, PI, 3/2 PI.
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
    ncols = L // l
    nrows = L // l
    x, y, phi = defect.create_defect_pair(d0, rho0, L, rot_angle=rot_angle)
    num, svx, svy = defect.coarse_grain(x, y, phi, L, L, l)
    VMpy.ini(x, y, np.cos(phi), np.sin(phi), seed, v0, eta, eps, L, L)
    FFMpegWriter = animation.writers["ffmpeg"]
    writer = FFMpegWriter(fps=2, metadata=dict(artist="Matplotlib"))
    filename = "data/%d_%g_%g_%g_%d.mp4" % (L, rho0, eta, eps, seed)
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))
    with writer.saving(fig, filename, dpi=100):
        im0 = axes[0].imshow(
            np.zeros((L, L)),
            animated=True,
            extent=[0, L, 0, L],
            origin="lower",
            vmin=0,
            vmax=rho0 * 2.5)
        im1 = axes[1].imshow(
            np.zeros((L, L)),
            animated=True,
            extent=[0, L, 0, L],
            origin="lower",
            vmin=-np.pi,
            vmax=np.pi,
            cmap="hsv")
        plt.colorbar(im0, ax=axes[0], orientation="horizontal")
        plt.colorbar(im1, ax=axes[1], orientation="horizontal")
        axes[0].set_title("density")
        axes[1].set_title("orientation")
        title = plt.suptitle("")
        title_template = r"$L=%d,\ \rho_0=%g, \eta=%g,\ \epsilon=%g,\ d_0=%d,\ " \
            + r"t=%d,\ \langle\phi\rangle=%.5f$"
        plt.tight_layout(rect=[0, 0, 1, 0.96])
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
            phi2d = np.arctan2(svy, svx).reshape(nrows, ncols)
            if t < 100:
                ymin, ymax = nrows * 3 // 8, nrows * 5 // 8
                xmin, xmax = ncols * 3 // 8, ncols * 5 // 8
            elif t < 500:
                ymin, ymax = nrows * 1 // 4, nrows * 3 // 4
                xmin, xmax = ncols * 1 // 4, ncols * 3 // 4
            else:
                ymin, ymax = 0, nrows
                xmin, xmax = 0, ncols
            im0.set_data(num.reshape(nrows, ncols)[ymin:ymax, xmin:xmax])
            im1.set_data(phi2d[ymin:ymax, xmin:xmax])
            title.set_text(title_template % (L, rho0, eta, eps, d0, t,
                                             order_para))
            im0.set_extent([ymin, ymax, xmin, xmax])
            im1.set_extent([ymin, ymax, xmin, xmax])
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
