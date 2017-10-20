import VMpy
import numpy as np
import matplotlib.pyplot as plt


if __name__ == "__main__":
    nBird = 10000
    np.random.seed(1)
    x, y, vx, vy = np.random.rand(4, nBird)
    Lx = 100
    Ly = 100
    x *= Lx
    y *= Ly
    vx *= 2 * np.pi
    vy *= 2 * np.pi

    VMpy.ini(x, y, vx, vy, 123, 0.5, 0.35, Lx, Ly)
    VMpy.run(10000)

    l = 1
    ncols = int(Lx / l)
    nrows = int(Ly / l)
    ncells = ncols * nrows
    num = np.zeros(ncells, int)
    vx = np.zeros(ncells)
    vy = np.zeros(ncells)

    VMpy.get_coarse_grained_snap(num, vx, vy, l)
    theta = np.arctan2(vy, vx)
    plt.imshow(theta.reshape(nrows, ncols), origin="lower", cmap="hsv")
    plt.colorbar()
    plt.show()
    plt.close()
