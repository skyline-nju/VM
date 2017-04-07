#ifndef IO_DATA_H
#define IO_DATA_H
#include <functional>
#include <vector>
#include <ctime>
#include "node.h"
#include "cmdline.h"

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


class CoarseGrain
{
public:
  CoarseGrain(double eta,
              double eps,
              double Lx,
              double Ly,
              int nBird,
              unsigned long long seed,
              const cmdline::parser &cmd,
              std::ofstream &log,
              bool &flag);
  void write(const Node *bird, int nBird, int step);
  void set_cell_size(double Lx, double Ly, const cmdline::parser &cmd);
  void set_output(double eta, double eps, double Lx, double Ly, int nBird,
                  unsigned long long seed, const cmdline::parser &cmd);
  void coarse_grain(const Node *bird, int nBird,
                    int *count, float *vx, float *vy,
                    double &svx, double &svy) const;
  void save_as_Bbb_format(
      const int *_count, const float *_vx, const float *_vy);
private:
  std::ofstream fout;
  std::vector<int> vec_frames;
  std::string format;
  char filename[100];
  int idx_cur_frame;
  int ncells;
  int ncols;
  int nrows;
  double lx;
  double ly;
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
};
#endif
