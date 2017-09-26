#ifndef IO_DATA_H
#define IO_DATA_H
#include <functional>
#include <vector>
#include <ctime>
#include "node.h"
#include "cmdline.h"
#include "spatial_corr.h"

class OrderPara
{
public:
  OrderPara(double eta,
            double epsilon,
            double rho0,
            double Lx,
            double Ly,
            unsigned long long seed,
            int phi_dt,
            std::ofstream &log);
  void out(const Node *bird, int nBird, int step);
private:
  std::ofstream fout;
  int dt;
};

class InSnapshot
{
public:
  InSnapshot(const std::string file);
  ~InSnapshot();
  void read();
  void read_block(int idx_frame);
  void from_file(int idx_frame, double &_Lx, double &_Ly,
                 std::vector<float> &x, std::vector<float> &y,
                 std::vector<float> &theta);
private:
  std::ifstream fin;
  bool is_read_block;
  double Lx;
  double Ly;
  int nBird;
  float *buff;
};


class OutSnapshot
{
public:
  OutSnapshot(double _eta,
              double _eps,
              double _rho0,
              double _Lx,
              double _Ly,
              unsigned long long _seed,
              int nBird,
              const cmdline::parser &cmd,
              std::ofstream &log,
              bool &flag);
  void write(const Node *bird, int nBird);
  void to_file(const Node *bird, int nBird, int step);
private:
  std::ofstream fout;
  bool is_one_file;
  int out_dt;
  double eta;
  double eps;
  double rho0;
  double Lx;
  double Ly;
  unsigned long long int seed;
};


// output coarse-grained snapshots
class CoarseGrain
{
public:
  CoarseGrain(int nBird, const cmdline::parser &cmd,
              std::ofstream &log, bool &flag);
  void write(const Node *bird, int nBird, int step);
  void set_box_size(const cmdline::parser &cmd);
  void set_output(const cmdline::parser &cmd, int nBird);

private:
  std::ofstream fout;
  std::vector<int> vec_frames;
  std::string format;
  char filename[100];
  int idx_cur_frame;
  size_t ncells;
  int ncols;
  int nrows;
  double lx;
  double ly;
};


// output correlation functions of density and velocity averaged spherially
class Corr_r
{
public:
  Corr_r(const cmdline::parser &cmd, int nBird);
  ~Corr_r();
  void instant(const Node *bird, int nBird);
  void output(const Node *bird, int nBird, int t);
private:
  std::ofstream fout;
  std::vector<int> frames;
  std::vector<int>::iterator iter;
  std::vector<double> r;
  SpatialCorr2d *corr2d;
  CircleAverage *circle_ave;
  double lBox;
  int ncols;
  size_t ncells;
  double *c_rho_r;
  double *c_v_r;
  double rho_m;
  double vx_m;
  double vy_m;
};


class Output
{
public:
  Output(double rho0, double Ly, int nBird, const cmdline::parser &cmd);
  ~Output();
  void out(const Node *bird, int nBird, int step);
  void ini_fout(double eta, double eps, double rho0,
                double Lx, double Ly, unsigned long long seed);
private:
  int interval;
  std::ofstream fout;
  std::time_t beg_time;
  std::vector<std::function<void(const Node*, int, int)>> fout_vec;
  OrderPara *phi;
  OutSnapshot *snap;
  CoarseGrain *cg;
  Corr_r *cr;
};

template <typename T1, typename T2>
void coarse_grain(const Node *bird, int nBird, T1 * num, T2 * vx, T2 * vy,
                  size_t ncells, int ncols, int nrows, double lx, double ly) {
  for (size_t i = 0; i < ncells; i++) {
    num[i] = 0;
    vx[i] = 0;
    vy[i] = 0;
  }
  for (unsigned int i = 0; i < nBird; i++) {
    int col = int(bird[i].x / lx);
    if (col < 0 || col >= ncols) col = 0;
    int row = int(bird[i].y / ly);
    if (row < 0 || row >= nrows) row = 0;
    int j = col + ncols * row;
    num[j]++;
    vx[j] += bird[i].vx;
    vy[j] += bird[i].vy;
  }
  for (size_t i = 0; i < ncells; i++) {
    if (num[i] > 0) {
      vx[i] /= num[i];
      vy[i] /= num[i];
    }
  }
}

template <typename T>
void coarse_grain(const Node *bird, int nBird, T *num,
                  size_t ncells, int ncols, int nrows, double lx, double ly,
                  bool flag_filter) {
  for (size_t i = 0; i < ncells; i++) {
    num[i] = 0;
  }
  for (unsigned i = 0; i < nBird; i++) {
    if (!flag_filter || bird[i].vx > 0) {
      int col = int(bird[i].x / lx);
      if (col < 0 || col >= ncols) col = 0;
      int row = int(bird[i].y / ly);
      if (row < 0 || row >= nrows) row = 0;
      num[col + ncols* row]++;
    }
  }
}
#endif
