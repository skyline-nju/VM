import sys
import numpy as np
from gsd import hoomd
from read_gsd import read_one_frame


def duplicate(s: hoomd.Snapshot, nx: int, ny: int) -> hoomd.Snapshot:
    N = s.particles.N * nx * ny
    lx = s.configuration.box[0]
    ly = s.configuration.box[1]
    Lx, Ly = lx * nx, ly * ny
    pos = np.zeros((N, 3), dtype=np.float32)
    for j in range(ny):
        for i in range(nx):
            beg = (j * nx + i) * s.particles.N
            end = beg + s.particles.N
            pos[beg:end, 0] = s.particles.position[:, 0] + lx / 2 + i * lx
            pos[beg:end, 1] = s.particles.position[:, 1] + ly / 2 + j * ly
            pos[beg:end, 2] = s.particles.position[:, 2]
    pos[:, 0] -= Lx / 2
    pos[:, 1] -= Ly / 2
    s2 = hoomd.Snapshot()
    s2.configuration.box = [Lx, Ly, 1, 0, 0, 0]
    s2.particles.N = N
    s2.particles.position = pos
    s2.configuration.step = 0
    return s2


def scale(s: hoomd.Snapshot, nx: int, ny: int, eps=0) -> hoomd.Snapshot:
    N = s.particles.N * nx * ny
    lx = s.configuration.box[0]
    ly = s.configuration.box[1]
    Lx, Ly = lx * nx, ly * ny
    pos = np.zeros((N, 3), dtype=np.float32)
    for i in range(nx * ny):
        beg = i * s.particles.N
        end = beg + s.particles.N
        pos[beg:end, 0] = s.particles.position[:, 0] * nx
        pos[beg:end, 1] = s.particles.position[:, 1] * ny
        pos[beg:end, 2] = s.particles.position[:, 2]
    if nx > 1:
        pos[:, 0] += (np.random.rand(N) - 0.5) * eps * nx
        mask = pos[:, 0] < Lx/2
        pos[:, 0][mask] += Lx
        mask = pos[:, 0] >= Lx/2
        pos[:, 0][mask] -= Lx
    # if ny > 1:
        pos[:, 1] += (np.random.rand(N) - 0.5) * eps * ny
        mask = pos[:, 1] < Ly/2
        pos[:, 1][mask] += Ly
        mask = pos[:, 1] >= Ly/2
        pos[:, 1][mask] -= Ly
    s2 = hoomd.Snapshot()
    s2.configuration.box = [Lx, Ly, 1, 0, 0, 0]
    s2.particles.N = N
    s2.particles.position = pos
    s2.configuration.step = 0
    return s2


def adjust_density(s: hoomd.Snapshot, phi: float, mode="copy") -> hoomd.Snapshot:
    def add_new_particles(n):
        pos_new = np.zeros((n, 3), dtype=np.float32)
        pos_new[:, 0] = (np.random.rand(n) - 0.5) * Lx
        pos_new[:, 1] = (np.random.rand(n) - 0.5) * Ly
        pos_new[:, 2] = (np.random.rand(n) - 0.5) * np.pi * 2
        return pos_new

    def extend_pos(pos_new, pos_old, mode):
        n_new = pos_new.shape[0]
        n_old = pos_old.shape[0]
        m = n_new // n_old
        n_res = n_new - m * n_old
        end = 0
    
        for i in range(m):
            beg = i * n_old
            end = beg + n_old
            pos_new[beg: end] = pos_old
        
        if n_res > 0:
            if mode == "rand":
                pos_new[end:] = add_new_particles(n_res)
            elif mode == "copy":
                np.random.shuffle(pos_old)
                pos_new[end:] = pos_old[:n_res]
            else:
                print("Error, mode should be copy or rand")
        

    Lx = int(s.configuration.box[0])
    Ly = int(s.configuration.box[1])

    n = int(round(Lx * Ly * phi))

    pos0 = s.particles.position
    n0 = s.particles.N
    pos = np.zeros((n, 3), dtype=np.float32)
    if n0 < n:
        extend_pos(pos, pos0, mode)
        print("add %d A particles" % (n - n0))
    elif n0 > n:
        np.random.shuffle(pos)
        pos = pos0[:n]
        print("remove %d A particles" % (n0 - n))
    
    s2 = hoomd.Snapshot()
    s2.configuration.box = [Lx, Ly, 1, 0, 0, 0]
    s2.particles.N = n
    s2.particles.position = pos
    s2.configuration.step = 0
    return s2
    

def inverse(x, y, theta, xc, yc):
    x_inv = 2 * xc - x
    y_inv = 2 * yc - y
    theta_inv = theta + np.pi
    theta_inv[theta_inv > np.pi] -= np.pi * 2
    theta_inv[theta_inv < -np.pi] += np.pi * 2
    return x_inv, y_inv, theta_inv


def shift_pos(s: hoomd.Snapshot, dx: float, dy: float) -> hoomd.Snapshot:
    Ly = s.configuration.box[1]
    Lx = s.configuration.box[0]

    x = s.particles.position[:, 0]
    y = s.particles.position[:, 1]

    x += dx
    y += dy
    if dy > 0:
        y[y >= Ly / 2] -= Ly
    else:
        y[y < -Ly / 2] += Ly
    if dx > 0:
        x[x >= Lx / 2] -= Lx
    else:
        x[x < -Lx / 2] += Lx
    s.particles.position[:, 0] = x
    s.particles.position[:, 1] = y
    return s


if __name__ == "__main__":
    folder = "AsymVM/data"
    basename = "L512_128_e0.03_r0.65_a0.1_s2000.gsd"

    fname = f"{folder}/{basename}"
    snap = read_one_frame(fname, -1)

    snap = adjust_density(snap, 3)
   
    snap.configuration.step = 0

    fout = f"{folder}/L512_128_e0.03_r3_a0.1_s2001.gsd"
    f = hoomd.open(name=fout, mode='wb')
    f.append(snap)



