import VMpy
import numpy as np
import matplotlib.pyplot as plt


if __name__ == "__main__":
    nBird = 100000
    np.random.seed(1)
    x, y, vx, vy = np.random.rand(4, nBird)
    Lx = 200
    Ly = 200
    x *= Lx
    y *= Ly
    vx *= 2 * np.pi
    vy *= 2 * np.pi

    VMpy.ini(x, y, vx, vy, 123, 0.5, 0.35, Lx, Ly)

    l = 1
    ncols = int(Lx / l)
    nrows = int(Ly / l)
    ncells = ncols * nrows
    num = np.zeros(ncells, int)
    vx = np.zeros(ncells)
    vy = np.zeros(ncells)

    for i in range(10):
        VMpy.run(1000)
        VMpy.get_coarse_grained_snap(num, vx, vy, l)
        plt.imshow(num.reshape(Lx, Ly), origin="lower", vmin=0, vmax=5)
        plt.colorbar()
        plt.show()
        plt.close()
