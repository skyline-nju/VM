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


def animate(a, L, eta, N, dN, psi=0, seed=1, dt=0.05):
    phi = create_defect_pair(a, L, psi)
    XY.ini(L, dt, eta, seed, phi)
    FFMpegWriter = animation.writers["ffmpeg"]
    writer = FFMpegWriter(fps=10, metadata=dict(artist="Matplotlib"))
    filename = "data/L=%d_Cl=%g.mp4" % (L, eta)
    fig = plt.figure()
    im = plt.imshow(
        np.zeros((L, L)),
        animated=True,
        extent=[0, L, 0, L],
        origin="lower",
        vmin=0,
        vmax=1,
        cmap="Greys")
    title = plt.title("", fontsize="xx-large")
    title_template = r"$L=%d,\ t=%g,$ order parameter: %.4f"
    cb = plt.colorbar()
    cb.set_label(r"$\sin ^2 (2 \phi)$")
    with writer.saving(fig, filename, dpi=100):
        step = 0
        while step < N:
            PHY = np.sin(2 * phi.reshape(L, L))**2
            im.set_data(PHY.reshape(L, L))
            order_para = np.sqrt(
                np.mean(np.cos(phi))**2 + np.mean(np.sin(phi))**2)
            title.set_text(title_template % (L, step * dt, order_para))
            writer.grab_frame()
            XY.run(dN, phi)
            print("N=%d" % step)
            step += dN


if __name__ == "__main__":
    animate(25, 512, 2, 30000, 20)
