#pragma once
#include <vector>
#include "vect.h"
#include "configure.h"

template <typename TPar>
void cal_order_para(const std::vector<TPar> &p_arr, double &phi, double &theta) {
  double svx = 0;
  double svy = 0;
  auto end = p_arr.cend();
  for (auto it = p_arr.cbegin(); it != end; ++it) {
    svx += (*it).vx;
    svy += (*it).vy;
  }
  double vx_m = svx / p_arr.size();
  double vy_m = svy / p_arr.size();
  phi = sqrt(vx_m * vx_m + vy_m * vy_m);
  theta = std::atan2(vy_m, vx_m);
}

void ini_order_para_exporter(double eta, double eps, double L, unsigned long long seed);

void output_order_para(double phi, double theta);

void ini_XY_exporter(double eta, double eps, double L, unsigned long long seed);

void output_XY_frame_head(const Vec_2<double> &l, int t, int n_par);

void output_XY_frame_tail();

void output_XY_one(double x, double y, double theta);

void ini_traj(double eta, double eps, double L, unsigned long long seed);

void output_traj(const std::vector<double> &x_arr, const std::vector<double> &y_arr, int t);

template <typename TPar>
void output_traj(const std::vector<TPar> &p_arr, int t, const Vec_2<double> &l) {
  int n = 4;
  std::vector<double> x_arr;
  std::vector<double> y_arr;
  x_arr.reserve(n);
  y_arr.reserve(n);
  for (int i = 0; i < n; i++) {
    x_arr.push_back(p_arr[i].x + l.x * p_arr[i].offset.x);
    y_arr.push_back(p_arr[i].y + l.y * p_arr[i].offset.y);
  }
  output_traj(x_arr, y_arr, t);
}

template <typename TPar>
void output_XY(const std::vector<TPar> &p_arr, int t, const Vec_2<double> &l) {
  output_XY_frame_head(l, t, p_arr.size());
  auto end = p_arr.cend();
  for (auto it = p_arr.cbegin(); it != end; ++it) {
    double x = (*it).x;
    double y = (*it).y;
    double theta = atan2((*it).vx, (*it).vy);
    output_XY_one(x, y, theta);
  }
  output_XY_frame_tail();
}

template <typename TPar>
void coarse_grain(const std::vector<TPar> &p_arr, unsigned short *num, float *vx, float *vy) {
  auto end = p_arr.cend();
  for (auto it = p_arr.cbegin(); it != end; ++it) {
    num[(*it).cell_index] += 1;
    vx[(*it).cell_index] += (*it).vx;
    vy[(*it).cell_index] += (*it).vy;
  }
}

void output_coarse_grained_snap(unsigned short *num, float *vx, float *vy, int n, int t,
                                double eta, double eps, double L, unsigned long long seed);

template <typename TPar>
void output_coarse_grained_snap(const std::vector<TPar> &p_arr, int t, const Vec_2<double> &l,
                                double eta, double eps, unsigned long long seed) {
  int n = int(l.x) * int(l.y);
  unsigned short *num = new unsigned short[n] {};
  float *vx = new float[n] {};
  float *vy = new float[n] {};
  coarse_grain(p_arr, num, vx, vy);
  output_coarse_grained_snap(num, vx, vy, n, t, eta, eps, l.x, seed);
  delete[] num;
  delete[] vx;
  delete[] vy;
}

class MSD {
public:
  template <typename TPar>
  MSD(const std::vector<TPar> &p_arr, int t, const Vec_2<double> &l, int sep);

  template <typename TPar>
  void cal_r_square(const std::vector<TPar> &p_arr, int t, const Vec_2<double> &l);

  template <typename TPar>
  int record(const std::vector<TPar> &p_arr, int t, const Vec_2<double> &l);

  void output(double eta, double eps, double L, unsigned long long seed);

  const std::vector<double> &r_square() const { return r_square_; }

  const std::vector<int> &t_arr() const { return t_arr_; }
protected:
  std::vector<Vec_2<double>> coord0;
  std::vector<double> r_square_;
  std::vector<int> t_arr_;
  int t0_;
  int sep_;
};

template <typename TPar>
MSD::MSD(const std::vector<TPar>& p_arr, int t, const Vec_2<double> &l, int sep) {
  coord0.reserve(p_arr.size());
  auto end = p_arr.cend();
  for (auto it = p_arr.cbegin(); it != p_arr.cend(); ++it) {
    double x = (*it).x + l.x * (*it).offset.x;
    double y = (*it).y + l.y * (*it).offset.y;
    coord0.emplace_back(x, y);
  }
  t0_ = t;
  sep_ = sep;
}

template <typename TPar>
void MSD::cal_r_square(const std::vector<TPar>& p_arr, int t, const Vec_2<double>& l) {
  double x2_sum = 0;
  double y2_sum = 0;
  size_t n = p_arr.size();
  for (int i = 0; i < n; i++) {
    const double dx = p_arr[i].x + p_arr[i].offset.x * l.x - coord0[i].x;
    const double dy = p_arr[i].y + p_arr[i].offset.y * l.y - coord0[i].y;
    x2_sum += dx * dx;
    y2_sum += dy * dy;
  }
  double r2 = (x2_sum + y2_sum) / n;
  r_square_.push_back(r2);
  t_arr_.push_back(t - t0_);
}

template <typename TPar>
int MSD::record(const std::vector<TPar>& p_arr, int t, const Vec_2<double>& l) {
  int flag;
  if ( t > t0_ && (t - t0_) % sep_ == 0) {
    cal_r_square(p_arr, t, l);
    flag = 1;
  } else {
    flag = 0;
  }
  return flag;
}

void output_mean_msd(const std::vector<MSD> &msd_arr,
                     double eta, double eps, double L, unsigned long long seed);

void ini_output_msd(double eta, double eps, double L, unsigned long long seed, int t0);

void output_msd(const MSD &msd, int idx);