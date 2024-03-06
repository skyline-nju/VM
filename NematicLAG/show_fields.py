import numpy as np
import matplotlib.pyplot as plt
import struct



if __name__ == "__main__":

    Lx = 256
    Ly = 128
    rho0 = 1.85
    seed = 7
    fin = "data/L%d_%d_b2_e0.9_r%g_D1_s%d_dt10000_t0.bin" % (Lx, Ly, rho0, seed)



    frame_size = Lx * Ly * 4 * 2
    with open(fin, "rb") as f:
        f.seek(0, 2)
        filesize = f.tell()

        n_frames = filesize//frame_size
        print("find", n_frames, "frames")
        f.seek(frame_size * (n_frames - 20))
        # f.seek(0)

        while f.tell() < filesize:
            buf = f.read(frame_size)
            data = np.array(struct.unpack("%dH" % (Lx * Ly * 4), buf))
            data = data.reshape(Lx * Ly, 4)
            n_p0, n_m0, n_0p, n_0m = data[:, 0].reshape(Ly, Lx), data[:, 1].reshape(Ly, Lx), data[:, 2].reshape(Ly, Lx), data[:, 3].reshape(Ly, Lx)
            rho = n_p0 + n_m0 + n_0p + n_0m
            Q_xx = n_p0 + n_m0
            Q_yy = n_0p + n_0m

            fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(9, 4), constrained_layout=True, sharex=True, sharey=True)
            im1 = ax1.imshow(rho, origin="lower", vmax=6)

            mask = rho > 0
            q_xx = Q_xx * 1.0
            q_yy = Q_yy * 1.0
            q_xx[mask] /= rho[mask]
            q_yy[mask] /= rho[mask]
            im2 = ax2.imshow(q_xx, origin="lower", vmin=0, vmax=1)
            im3 = ax3.imshow(q_yy, origin="lower", vmin=0, vmax=1)

            cb1 = fig.colorbar(im1, ax=ax1, orientation="horizontal")
            cb2 = fig.colorbar(im2, ax=ax2, orientation="horizontal")
            cb3 = fig.colorbar(im3, ax=ax3, orientation="horizontal")
            cb1.set_label(r"$\rho$")
            cb2.set_label(r"$Q_{xx}/\rho$")
            cb3.set_label(r"$Q_{yy}/\rho$")

            plt.show()
            plt.close()
            

