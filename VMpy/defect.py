'''
    Study the evolusion of defect pairs in the Vicsek model.
'''

import numpy as np
import matplotlib.pyplot as plt
import VMpy


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
    for k in range(x.size):
        i = int(x[k])
        j = int(y[k])
        vx_new[j, i] += vx[k]
        vy_new[j, i] += vy[k]
    phi_new = np.arctan2(vy_new, vx_new)
    return phi_new


if __name__ == "__main__":
    d0 = 20
    rho0 = 4
    L = 200
    seed = 123
    v0 = 0
    eta = 0.05
    l = 1
    ncols = L // l
    nrows = L // l
    num = np.zeros(ncols * nrows, int)
    svx = np.zeros(ncols * nrows)
    svy = np.zeros(ncols * nrows)

    x, y, phi = create_defect_pair(d0, rho0, L, rot_angle=-0.5 * np.pi)
    phi2d = coarse_grain(x, y, phi, L, L, l)
    plt.imshow(phi2d, origin="lower", cmap="hsv", interpolation="none")
    plt.title(r"$t=%d$" % 0)
    plt.colorbar()
    # plt.show()
    plt.savefig("data/%g_%d_%d_%d.png" % (eta, v0, L, 0))
    plt.close()

    VMpy.ini(x, y, np.cos(phi), np.sin(phi), seed, v0, eta, L, L)

    dt = 1000
    for i in range(100):

        VMpy.run(dt)
        VMpy.get_coarse_grained_snap(num, svx, svy, l)
        order_para = np.sqrt(np.sum(svx)**2 + np.sum(svy)**2) / (L**2 * rho0)
        phi2d = np.arctan2(svy, svx).reshape(nrows, ncols)

        plt.imshow(phi2d, origin="lower", cmap="hsv", interpolation="none")
        plt.title(r"$t=%d\, \langle\phi\rangle=%.4f$" % ((i + 1) * dt,
                                                         order_para))
        plt.colorbar()
        # plt.show()
        plt.savefig("data/%g_%d_%d_%d.png" % (eta, v0, L, (i+1) * dt))
        plt.close()
