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
    import matploblib.animation as animation
    dest = "/ffmpeg-3.3-64bit-static/ffmpeg"
    path1 = "/home-gk/users/nscc1185/Applications"
    path2 = "/home-yw/users/nsyw449_YK/dy/Program"
    if os.path.exist(path1):
        plt.rcParams['animation.ffmpeg_path'] = path1 + dest
    elif os.path.exist(path2):
        plt.rcParams['animation.ffmpeg_path'] = path2 + dest
    else:
        print("Error, cannot find ffmpeg")
        sys.exit()


def make_movie(L, rho0, d0, seed, v0, eta, rot_angle, nstep, interval):
    l = 1
    ncols = L // l
    nrows = L // l
    x, y, phi = defect.create_defect_pair(d0, rho0, L, rot_angle=rot_angle)
    num, svx, svy = defect.coarse_grain(x, y, phi, L, L, l)
    VMpy.ini(x, y, np.cos(phi), np.sin(phi), seed, v0, eta, L, L)
    FFMpegWriter = animation.writers["ffmpeg"]
    writer = FFMpegWriter(fps=2, metadata=dict(artist="Matplotlib"))
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))
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
        vmax=np.pi)
    plt.colorbar(im0, ax=axes[0], orientation="horizontal")
    plt.colorbar(im1, ax=axes[1], orientation="horizontal")
    title = plt.suptitle("")
    title_template = r"$t=%d,\ \langle\phi\rangle=%.4f$"
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    filename = "data/%d_%g_%g_%g_%d.mp4" % (L, rho0, eta, v0, seed)
    with writer.saving(fig, filename, dpi=100):
        t = 0
        while t <= nstep:
            print("t=%d" % t)
            if t > 0:
                VMpy.run(interval)
                VMpy.get_coarse_grained_snap(num, svx, svy, l)
            t += interval
            order_para = np.sqrt(
                np.sum(svx)**2 + np.sum(svy)**2) / (L**2 * rho0)
            phi2d = np.arctan2(svy, svx).reshape(nrows, ncols)
            im0.set_data(num.reshape(nrows, ncols))
            im1.set_data(phi2d)
            title.set_text(title_template % (t, order_para))
            writer.grab_frame()


if __name__ == "__main__":
    import defect
    import VMpy
    d0 = 10
    rho0 = 4
    L = 1024
    seed = 124
    v0 = 0.5
    eta = 0.2
    l = 1
    make_movie(L, rho0, d0, seed, v0, eta, 0, 2000, 10)
