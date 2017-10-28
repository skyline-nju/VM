""" Test XYmodel."""

import numpy as np
import XYmodel as XY
import sys
import os
import platform
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


def create_defect_pair(a, L, rot_angle=0):
    """ Create a pair of defects with separation 2a.

    Parameters:
    --------
    a : int
        Half of separation between two defects.
    L : int
        System size.
    rot_angle: double, optional
        Rotate the global orientational fields by this angle.
    """
    x = np.arange(L) - 0.5 * L + 0.5
    y = np.arange(L) - 0.5 * L + 0.5
    xx, yy = np.meshgrid(x, y)
    phi = np.arctan2(yy, xx - a) - np.arctan2(yy, xx + a) + rot_angle
    phi[phi > np.pi] -= np.pi * 2
    phi[phi < -np.pi] += np.pi * 2
    return phi.flatten()


def create_rand_uniaxial_field(L):
    mask = np.random.rand(L * L) < 0.5
    psi = np.zeros(L * L)
    psi[mask] = np.pi
    return psi


def animate(L, eta, N, dN, phi=None, h=None, seed=1, dt=0.05, cmap="Greys"):
    if phi is None:
        phi = np.random.rand(L * L) * np.pi * 2
    if h is None:
        XY.ini(L, dt, eta, seed, phi)
        filename = "data/L=%d_Cl=%g.mp4" % (L, eta)
        title_template = r"$L=%d,\ C_L=%g,\ t=%g,$ order parameter: %.4f"
    else:
        psi = create_rand_uniaxial_field(L)
        XY.ini(L, dt, eta, seed, phi, h, psi)
        filename = "data/L=%d_Cl=%g_h=%g.mp4" % (L, eta, h)
        title_template = r"$L=%d,\ C_L=%g,\ h=%g,\ t=%g,$" \
            + "\norder parameter: %.4f"
    FFMpegWriter = animation.writers["ffmpeg"]
    writer = FFMpegWriter(fps=10, metadata=dict(artist="Matplotlib"))
    fig = plt.figure()
    if cmap == "Greys":
        vmin = 0
        vmax = 1
    else:
        vmin = -np.pi
        vmax = np.pi
    im = plt.imshow(
        np.zeros((L, L)),
        animated=True,
        extent=[0, L, 0, L],
        origin="lower",
        vmin=vmin,
        vmax=vmax,
        cmap=cmap)
    title = plt.title("", fontsize="xx-large")
    cb = plt.colorbar()
    cb.set_label(r"$\sin ^2 (2 \phi)$")
    with writer.saving(fig, filename, dpi=100):
        step = 0
        while step < N:
            order_para = np.sqrt(
                np.mean(np.cos(phi))**2 + np.mean(np.sin(phi))**2)
            if h is None:
                PHY = np.sin(2 * phi.reshape(L, L))**2
                im.set_data(PHY.reshape(L, L))
                Title = title_template % (L, eta, step * dt, order_para)
            else:
                phi[phi < -np.pi] += 2 * np.pi
                phi[phi > np.pi] -= 2 * np.pi
                im.set_data(phi.reshape(L, L))
                Title = title_template % (L, eta, h, step * dt, order_para)
            title.set_text(Title)
            writer.grab_frame()
            XY.run(dN, phi)
            print("N=%d" % step)
            step += dN


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-L", type=int, help="System size")
    parser.add_argument("--Cl", type=float, help="Strength of noise")
    parser.add_argument(
        "-H", type=float, default=None, help="Strength of random field")
    parser.add_argument("--nStep", type=int, help="total steps to run")
    parser.add_argument("--dn", type=int, help="interval between two frames")
    parser.add_argument(
        "-d", type=float, default=None, help="Distance between two defects")
    parser.add_argument(
        "--spin_num", type=int, default=0, help="rotate multiple 90 degree")
    parser.add_argument("-s", type=int, default=0, help="random number seed")
    arg = parser.parse_args()
    if arg.s == 0:
        import datetime
        seed = int(datetime.datetime.now().strftime("%m%d%H%M%S"))
    else:
        seed = arg.s
    if arg.H is None:
        phi = create_defect_pair(arg.d, arg.L, arg.spin_num * np.pi / 2)
        animate(arg.L, arg.Cl, arg.nStep, arg.dn, phi=phi, seed=seed)
    else:
        animate(
            arg.L, arg.Cl, arg.nStep, arg.dn, h=arg.H, cmap="hsv", seed=seed)
