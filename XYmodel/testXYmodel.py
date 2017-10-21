""" Test XYmodel."""

import numpy as np
import matplotlib.pyplot as plt
import XYmodel as XY


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


if __name__ == "__main__":
    L = 512
    dt = 0.05
    eta = 1
    n = L * L
    seed = 123
    a = 50
    phi = create_defect_pair(a, L, 0)
    # plt.imshow(phi.reshape(L, L), cmap="hsv", origin="lower")
    psi = np.sin(2 * phi.reshape(L, L)) ** 2
    plt.imshow(psi, cmap="Greys", origin="lower")
    order_para = np.sqrt(np.mean(np.cos(phi))**2 + np.mean(np.sin(phi))**2)
    plt.title(r"$t=%d\, \langle \phi\rangle = %.4f$" % (0, order_para))
    # plt.show()
    plt.savefig("data/%d_%d_%d.png" % (L, eta, 0))
    plt.close()

    XY.ini(L, dt, eta, seed, phi)

    step = 1000
    for i in range(200):
        XY.run(step, phi)
        # plt.imshow(phi.reshape(L, L), cmap="hsv", origin="lower")
        psi = np.sin(2 * phi.reshape(L, L)) ** 2
        plt.imshow(psi, cmap="Greys", origin="lower")
        order_para = np.sqrt(np.mean(np.cos(phi))**2 + np.mean(np.sin(phi))**2)
        plt.title(r"$t=%d\, \langle \phi\rangle = %.4f$" % ((i + 1) * step,
                                                            order_para))
        # plt.show()
        plt.savefig("data/g%d_%d_%d.png" % (L, eta, (i + 1) * step))
        plt.close()
