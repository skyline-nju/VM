#ifndef SPATIAL_CORR_H
#define SPATIAL_CORR_H
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fftw3.h>

class AutoCorr2d
{
public:
  AutoCorr2d(int _nrow, int _ncol);
  ~AutoCorr2d();
  void autocorr2(double * in_new, double * out_new);
  void autocorr2(double * in_new, double * out_new, double *Sk);
private:
  int ntot;
  int ntot_complex;

  double * in;
  fftw_complex *inter;
  double * out;

  fftw_plan fft;
  fftw_plan rfft;
};

class SpatialCorr2d: public AutoCorr2d
{
public:
  SpatialCorr2d(int ncols0, int nrows0, double Lx, double Ly);
  void eval(const int *num, const double *vx, const double *vy,
            double *c_rho, double *c_v,
            double &rho_m, double &vx_m, double &vy_m);
  void eval(const int *num, const double *vx, const double *vy,
            double *c_rho, double *c_v, double *s_rho, double *s_v,
            double &rho_m, double &vx_m, double &vy_m);
private:
  int ncols;
  int nrows;
  size_t ncells;
  double lx;
  double ly;
  double dA;
};

class CircleAverage
{
public:
  CircleAverage(int _ncols, int _nrows, double l, std::vector<double> &r);
  void eval(const double *corr2d, double *corr_r) const;
private:
  int ncols;
  int nrows;
  int half_nrows;
  int half_ncols;
  unsigned int RR_max;
  int R_thresh;
  unsigned int RR_thresh;
  int idx_thresh;
  int size_r;
  std::vector<unsigned int> count;
};

#endif