#ifndef IO_DATA_H
#define IO_DATA_H
#include <fstream>
#include <ctime>
#include <chrono>
#include "comn.h"
#include "cmdline.h"
#include "vm.h"

void ini_output(const cmdline::parser &cmd, const VM *birds);

void output(int i, const VM* birds);

class _Writer {
public:
  _Writer(const cmdline::parser &cmd) : idx_frame(0) {}
  virtual void set_frames(const cmdline::parser &cmd) = 0;
  virtual void write(int i, const VM *bird, std::ofstream &fout) = 0;

protected:
  std::vector<int> frames;
  int idx_frame;
};

class OrderParaWriter : public _Writer {
public:
  OrderParaWriter(const cmdline::parser &cmd, std::ofstream &fout);
  void set_frames(const cmdline::parser &cmd);
  void write(int i, const VM* bird, std::ofstream &fout);
};

class LogWriter : public _Writer {
public:
  LogWriter(const cmdline::parser &cmd, std::ofstream &fout);
  void set_frames(const cmdline::parser &cmd);
  void write(int i, const VM *bird, std::ofstream &fout);
private:
  std::chrono::time_point<std::chrono::system_clock> t_start;
};

class CoarseGrainSnapWriter : public _Writer {
public:
  CoarseGrainSnapWriter(const cmdline::parser &cmd, std::ofstream &fout);
  void set_frames(const cmdline::parser &cmd);
  void write(int i, const VM* bird, std::ofstream &fout);

private:
  double l;
  int ncols;
  int ncells;
  double dt;
  bool is_continous;
};

template <typename T1, typename T2>
void coarse_grain(T1 *num_cg, T2 *vx_cg, T2 *vy_cg, double &vx_m, double &vy_m,
                  double l, int ncols, int ntot, const VM *birds) {
  for (int i = 0; i < ntot; i++) {
    num_cg[i] = 0;
    vx_cg[i] = 0;
    vy_cg[i] = 0;
  }
  vx_m = 0;
  vy_m = 0;
  int N = birds->get_num_birds();
  for (int i = 0; i < N; i++) {
    double x, y, vx, vy;
    birds->get_x(i, x, y);
    birds->get_v(i, vx, vy);
    int col = int(x / l);
    int row = int(y / l);
    int j = col + ncols * row;
    num_cg[j]++;
    vx_cg[j] += vx;
    vy_cg[j] += vy;
    vx_m += vx;
    vy_m += vy;
  }
  vx_m /= N;
  vy_m /= N;
}
#endif
