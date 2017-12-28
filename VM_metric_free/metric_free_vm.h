#ifndef METRIC_FREE_VM_H
#define METRIC_FREE_VM_H
#include <fstream>
#include <ctime>
#include <chrono>
#include <string>

#include "bird.h"
#include "comn.h"
#include "cmdline.h"

void ini(cmdline::parser &cmd);

void run(int n);

void end();

template <typename T1, typename T2>
void coarse_grain(T1 *num, T2 *vx, T2 *vy, double *v_mean, double l,
                  std::vector<Pair_P_I> &pos_idx_pair,
                  std::vector<Bird_wo_pos> &birds) {
  int nBird = pos_idx_pair.size();
  int ncols = int(Bird_wo_pos::Lx / l);
  int nrows = int(Bird_wo_pos::Ly / l);
  int ncells = ncols * nrows;
  v_mean[0] = 0;
  v_mean[1] = 0;
  for (int i = 0; i < ncells; i++) {
    num[i] = 0;
    vx[i] = 0;
    vy[i] = 0;
  }
  for (int i = 0; i < nBird; i++) {
    int col = int(pos_idx_pair[i].first.x() / l);
    int row = int(pos_idx_pair[i].first.y() / l);
    if (col < 0) {
      col += ncols;
    } else if (col >= ncols) {
      col -= ncols;
    }
    if (row < 0) {
      row += nrows;
    } else if (row >= nrows) {
      row -= nrows;
    }
    int idx = col + row * ncols;
    num[idx] += 1;
    vx[idx] += birds[i].vx;
    vy[idx] += birds[i].vy;
    v_mean[0] += birds[i].vx;
    v_mean[1] += birds[i].vy;
  }
  v_mean[0] /= nBird;
  v_mean[1] /= nBird;
}

template <typename T1, typename T2>
void output_coarse_grained_snap(int t, T1 *num, T2 *vx, T2 *vy,
                                int ncells, double l,
                                std::ofstream &fout,
                                std::vector<Pair_P_I> &pos_idx_pair,
                                std::vector<Bird_wo_pos> &birds) {
  fout.write((char *)&t, sizeof(int));
  double v_mean[2];
  coarse_grain(num, vx, vy, v_mean, l, pos_idx_pair, birds);
  fout.write((char *)v_mean, sizeof(double) * 2);
  fout.write((char *)num, sizeof(T1) * ncells);
  fout.write((char *)vx, sizeof(T2) * ncells);
  fout.write((char *)vy, sizeof(T2) * ncells);
}
#endif