import VMpy
import numpy as np
import platform
import matplotlib
if platform.system() is "Windows":
    import matplotlib.pyplot as plt
else:
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt


if __name__ == "__main__":
    nBird = 1000
    np.random.seed(1)
    x, y, vx, vy = np.random.rand(4, nBird)
    Lx = 100
    Ly = 100
    x *= Lx
    y *= Ly
    vx *= 2 * np.pi
    vy *= 2 * np.pi

    VMpy.ini(x, y, vx, vy, 123, 0.5, 0.35, Lx, Ly)

    l = 1.0
    ncols = int(Lx / l)
    nrows = int(Ly / l)
    ncells = ncols * nrows
    
    num = np.zeros(ncells, np.int32)  # np.int32 is necessary on Linux platform
    vx = np.zeros(ncells)
    vy = np.zeros(ncells)

    VMpy.run(10)
    VMpy.get_coarse_grained_snap(num, vx, vy, l)
    plt.imshow(num.reshape(nrows, ncols), origin="lower", vmin=0, vmax=5)
    plt.colorbar()
    if platform.system() is "Windows":
        plt.show()
    else:
        plt.savefig("test.png")
    plt.close()
