import VMpy
import numpy as np
import matplotlib.pyplot as plt


if __name__ == "__main__":
    nBird = 100
    x, y, vx, vy = np.random.rand(4, nBird)
    Lx = 100
    Ly = 100
    x *= Lx
    y *= Ly
    vx *= 2 * np.pi
    vy *= 2 * np.pi
    plt.plot(x, y, "o")
    plt.show()
    plt.close()
    VMpy.set_random_seed(123)
    VMpy.set_eta(0)
    VMpy.set_v0(0.5)
    VMpy.setLx(Lx)
    VMpy.setLy(Ly)
    VMpy.run(10, x, y, vx, vy)
    plt.plot(x, y, "o")
    plt.show()
    plt.close()
    VMpy.run(10, x, y, vx, vy)
    plt.plot(x, y, "o")
    plt.show()
    plt.close()
