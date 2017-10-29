import numpy as np
import cVMpy as cvm
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


def fly(nBird, eta, Lx, Ly, seed=124, dt=0.05, interval=20):
    def init():
        ax.set_xlim(0, Lx)
        ax.set_ylim(0, Ly)
        # title = ax.set_title("a")
        return title, ln

    def update(frame, *args):
        cvm.run(interval)
        xdata, ydata, theta = np.zeros((3, nBird))
        cvm.get_snap(xdata, ydata, theta)
        ln.set_data(xdata, ydata)
        # plt.title(r"$t=%g$" % (dt * frame))
        vxm = np.mean(np.cos(theta))
        vym = np.mean(np.sin(theta))
        phi = np.sqrt(vxm**2 + vym**2)
        title.set_text(r"$t=%g, \phi=%.3f, v_y / v_x=%.3f$" %
                       (dt * frame, phi, vym / vxm))
        return title, ln

    fig, ax = plt.subplots()
    ln, = plt.plot([], [], 'ro', animated=True, ms=2)
    title = ax.text(
        0.1,
        0.9,
        "",
        animated=True,
        transform=ax.transAxes,
        fontsize="x-large")
    # cvm.ini(xdata, ydata, theta, nBird, Lx, Ly, seed, dt)
    cvm.ini_rand(nBird, eta, Lx, Ly, seed, dt)
    cvm.set_h(0.2)
    cvm.set_period(4)

    frames = (np.arange(500) + 1) * interval
    ani = FuncAnimation(
        fig, update, init_func=init, frames=frames, repeat=False, blit=True)
    plt.show()
    plt.close()


if __name__ == "__main__":
    fly(8000, 0.2, 70, 70, dt=0.1, interval=40)
