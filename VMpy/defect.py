'''
    Study the evolusion of defect pairs in the Vicsek model.
'''

import numpy as np
import matplotlib.pyplot as plt


def show_defects(Lx, Ly, a, psi=0, ax=None, xc=0, yc=0, show_axis=False):
    if ax is None:
        ax = plt.subplot(111)
        flag_show = True
    else:
        flag_show = False
    x = np.arange(Lx) + 0.5 - 0.5 * Lx + xc
    y = np.arange(Ly) + 0.5 - 0.5 * Ly + yc
    xx, yy = np.meshgrid(x, y)
    phi1 = np.arctan2(yy - yc, xx - a - xc)
    phi2 = np.arctan2(yy - yc, xx + a - xc)
    phi = phi1 - phi2 + psi
    phi[phi < -np.pi] += 2 * np.pi
    phi[phi > np.pi] -= 2 * np.pi
    dx = np.cos(phi) * 0.5
    dy = np.sin(phi) * 0.5
    box = [x[0], x[-1], y[0], y[-1]]
    ax.imshow(phi, origin="lower", cmap="hsv", extent=box, alpha=0.7)
    ax.quiver(xx, yy, dx, dy, units='inches')
    ax.axis("scaled")
    if not show_axis:
        ax.axis("off")
    if flag_show:
        plt.show()
        plt.close()


def show_four_defect_pairs():
    Lx = 40
    Ly = 30
    a = 10
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12, 9))
    psi = [0, np.pi * 0.5, np.pi, 1.5 * np.pi]
    label = [
        r"(a) $\psi=0$", r"(b) $\psi=\pi / 2$", r"(c) $\psi= \pi$",
        r"(d) $\psi=3\pi/2$"
    ]
    for i, ax in enumerate(axes.flat):
        show_defects(Lx, Ly, a, psi[i], ax)
        ax.set_title(label[i], fontsize="xx-large")
    plt.tight_layout()
    # plt.show()
    plt.savefig("data/defect_pairs.pdf")
    plt.close()


def create_defect_pair(a, rho, Lx, Ly=None, rot_angle=0):
    """ Create a pair of defects with separation 2a.

    Parameters:
    --------
    a : int
        Half ot separation between two defects.
    rho : double
        Density of particles.
    Lx : int
        System size in the x direction.
    Ly : int, optional
        System size in the y direction.
    rot_angle ： double, optional
        Rotate the global orientational fields by this angle.
    """
    if Ly is None:
        Ly = Lx
    N = int(rho * Lx * Ly)
    x0, y0 = 0.5 * Lx, 0.5 * Ly
    x, y = np.random.rand(2, N)
    x *= Lx
    y *= Ly
    X = x - x0
    Y = y - y0
    phi1 = np.arctan2(Y, X - a)
    phi2 = np.arctan2(Y, X + a)
    phi = phi1 - phi2 + rot_angle
    phi[phi > np.pi] -= np.pi * 2
    phi[phi < -np.pi] += np.pi * 2
    return x, y, phi


def coarse_grain(x, y, phi, Lx, Ly, l=1):
    vx = np.cos(phi)
    vy = np.sin(phi)
    nx = Lx // l
    ny = Ly // l
    vx_new = np.zeros((ny, nx))
    vy_new = np.zeros((ny, nx))
    num = np.zeros((ny, nx), np.int32)  # np.int32 for Linux platform
    for k in range(x.size):
        i = int(x[k] / l)
        j = int(y[k] / l)
        vx_new[j, i] += vx[k]
        vy_new[j, i] += vy[k]
        num[j, i] += 1
    return num.flatten(), vx_new.flatten(), vy_new.flatten()


def plot_one_frame(num, svx, svy, t, axes=None):
    order_para = np.sqrt(np.sum(svx)**2 + np.sum(svy)**2) / (L**2 * rho0)
    phi2d = np.arctan2(svy, svx).reshape(nrows, ncols)

    if axes is None:
        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))
    im0 = axes[0].imshow(num.reshape(L, L), origin="lower", vmin=0, vmax=10)
    im1 = axes[1].imshow(
        phi2d,
        origin="lower",
        cmap="hsv",
        interpolation="none",
        vmin=-np.pi,
        vmax=np.pi)
    plt.colorbar(im0, ax=axes[0], orientation="horizontal")
    plt.colorbar(im1, ax=axes[1], orientation="horizontal")
    fig.suptitle(
        r"$t=%d\, \langle\phi\rangle=%.4f$" % (t, order_para),
        fontsize="x-large")
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    # plt.show()
    plt.savefig("data/%g_%d_%d_%d_%d.png" % (eta, v0, L, seed, t))
    plt.close()


if __name__ == "__main__":
    # import VMpy
    d0 = 10
    rho0 = 4
    L = 1024
    seed = 124
    v0 = 0.5
    eta = 0
    l = 1
    ncols = L // l
    nrows = L // l
    # num = np.zeros(ncols * nrows, int)
    # svx = np.zeros(ncols * nrows)
    # svy = np.zeros(ncols * nrows)
    # x, y, phi = create_defect_pair(d0, rho0, L, rot_angle=0)
    # num, svx, svy = coarse_grain(x, y, phi, L, L, l)
    # plot_one_frame(num, svx, svy, 0)

    # VMpy.ini(x, y, np.cos(phi), np.sin(phi), seed, v0, eta, L, L)

    # dt = 5
    # for i in range(200):

    #     VMpy.run(dt)
    #     VMpy.get_coarse_grained_snap(num, svx, svy, l)
    #     plot_one_frame(num, svx, svy, (i + 1) * dt)
    show_four_defect_pairs()
