import numpy as np
import matplotlib.pyplot as plt
import struct

def show_fields():
    Lx = 64
    Ly = 256
    beta = 4
    rho0 = 6
    eps = 0.9
    D = 0.1
    seed = 3001
    dt = 5000
    t_beg = 0
    fin = "data/L%d_%d_b%g_r%g_e%g_D%g_s%d_dt%d_t%d.bin" % (Lx, Ly, beta, rho0, eps, D, seed, dt, t_beg)


    frame_size = Lx * Ly * 4
    with open(fin, "rb") as f:
        f.seek(0, 2)
        filesize = f.tell()

        n_frames = filesize//frame_size
        print("find", n_frames, "frames")
        if n_frames >= 10:
            f.seek(frame_size * (n_frames - 10))
        else:
            f.seek(0)

        while f.tell() < filesize:
            buf = f.read(frame_size)
            sigma = np.array(struct.unpack("%dB" % (Lx * Ly * 4), buf)).reshape(Ly, Lx, 4).astype(int)
            rho = np.sum(sigma, axis=2)
            mx = sigma[:, :, 0] - sigma[:, :, 2]
            my = sigma[:, :, 1] - sigma[:, :, 3]
            # print(rho.shape, mx.shape)

            fig, axes = plt.subplots(1, 3, figsize=(9, 4), constrained_layout=True, sharex="col")
            im1 = axes[0].imshow(rho, origin="lower", vmin=0, vmax=10)

            mask = rho > 0

            px = mx.astype(float)
            px[mask] /= rho[mask]
            py = my.astype(float)
            py[mask] /= rho[mask]

            im2 = axes[1].imshow(px, origin="lower", vmin=-1, vmax=1, cmap="bwr")
            im3 = axes[2].imshow(py, origin="lower", vmin=-1, vmax=1, cmap="bwr")

            cb1 = fig.colorbar(im1, ax=axes[0], orientation="horizontal")
            cb2 = fig.colorbar(im2, ax=axes[1], orientation="horizontal")
            cb3 = fig.colorbar(im3, ax=axes[2], orientation="horizontal")

            cb1.set_label(r"$\rho$")
            cb2.set_label(r"$m_x/\rho$")
            cb3.set_label(r"$m_y/\rho$")

            plt.show()
            plt.close()


if __name__ == "__main__":
    Lx = 128
    Ly = 128
    beta = 3.2
    eta = -0.8
    rho0 = 10
    rho_thresh = 10
    eps = 1
    D = 0.1
    seed = 3001
    dt = 5000
    t_beg = 0
    fin = "data/L%d_%d_b%g_e%g_r%g_%g_v%g_D%g_s%d_dt%d_t%d.bin" % (Lx, Ly, beta, eta, rho0, rho_thresh, eps, D, seed, dt, t_beg)


    frame_size = Lx * Ly * 4
    with open(fin, "rb") as f:
        f.seek(0, 2)
        filesize = f.tell()

        n_frames = filesize//frame_size
        print("find", n_frames, "frames")
        if n_frames >= 10:
            f.seek(frame_size * (n_frames - 10))
        else:
            f.seek(0)

        while f.tell() < filesize:
            buf = f.read(frame_size)
            sigma = np.array(struct.unpack("%dB" % (Lx * Ly * 4), buf)).reshape(Ly, Lx, 4).astype(int)
            rho = np.sum(sigma, axis=2)
            mx = sigma[:, :, 0] - sigma[:, :, 2]
            my = sigma[:, :, 1] - sigma[:, :, 3]
            # print(rho.shape, mx.shape)

            fig, axes = plt.subplots(1, 3, figsize=(9, 4), constrained_layout=True, sharex="col")
            im1 = axes[0].imshow(rho, origin="lower", vmin=0, vmax=20)

            mask = rho > 0

            px = mx.astype(float)
            px[mask] /= rho[mask]
            py = my.astype(float)
            py[mask] /= rho[mask]

            im2 = axes[1].imshow(px, origin="lower", vmin=-1, vmax=1, cmap="bwr")
            im3 = axes[2].imshow(py, origin="lower", vmin=-1, vmax=1, cmap="bwr")

            cb1 = fig.colorbar(im1, ax=axes[0], orientation="horizontal")
            cb2 = fig.colorbar(im2, ax=axes[1], orientation="horizontal")
            cb3 = fig.colorbar(im3, ax=axes[2], orientation="horizontal")

            cb1.set_label(r"$\rho$")
            cb2.set_label(r"$m_x/\rho$")
            cb3.set_label(r"$m_y/\rho$")

            plt.show()
            plt.close()
            
