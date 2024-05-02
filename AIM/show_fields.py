import numpy as np
import matplotlib.pyplot as plt
import struct


if __name__ == "__main__":
    Lx = 512
    Ly = 32
    beta = 1.8
    rho0 = 3
    alpha = 0
    eps = 0.9
    D = 1
    seed = 4002
    dt = 1000
    t_beg = 0
    fin = "data/L%d_%d_b%g_r%g_a%g_e%g_D%g_s%d_dt%d_t%d.bin" % (Lx, Ly, beta, rho0, alpha, eps, D, seed, dt, t_beg)


    frame_size = Lx * Ly * 2 * 2
    with open(fin, "rb") as f:
        f.seek(0, 2)
        filesize = f.tell()

        n_frames = filesize//frame_size
        print("find", n_frames, "frames")
        f.seek(frame_size * (n_frames - 10))
        # f.seek(0)

        while f.tell() < filesize:
            buf = f.read(frame_size)
            data = np.array(struct.unpack("%dH%dh" % (Lx * Ly, Lx * Ly), buf)).reshape(2, Ly, Lx)
            rho, m = data

            fig, axes = plt.subplots(2, 2, figsize=(9, 4), constrained_layout=True, sharex="col")
            im1 = axes[1, 0].imshow(rho, origin="lower", vmin=0, vmax=5)

            mask = rho > 0
            phi = np.zeros_like(m)
            phi[mask] = m[mask]/rho[mask]

            im2 = axes[1, 1].imshow(phi, origin="lower", vmin=-1, vmax=1, cmap="bwr")

            cb1 = fig.colorbar(im1, ax=axes[1, 0], orientation="horizontal")
            cb2 = fig.colorbar(im2, ax=axes[1, 1], orientation="horizontal")
            cb1.set_label(r"$\rho$")
            cb2.set_label(r"$m/\rho$")

 

            m_x = np.mean(m, axis=0)
            x = np.arange(Lx) + 0.5
            rho_x = np.mean(rho, axis=0)
            axes[0, 0].plot(x, rho_x)
            axes[0, 1].plot(x, m_x)
            axes[0, 1].axhline(0, linestyle="dashed", c="k")

            plt.show()
            plt.close()
            
