#include "spatial_corr.h"
using namespace std;

AutoCorr2d::AutoCorr2d(int _nrow, int _ncol) {
  // shape of input/output array = (_nrows, _ncols)
  int ncol_complex = _ncol / 2 + 1;
  ntot = _nrow * _ncol;
  ntot_complex = _nrow * ncol_complex;
  in = (double *)fftw_malloc(sizeof(double) * ntot);
  inter = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * ntot_complex);
  out = (double *)fftw_malloc(sizeof(double) * ntot);
  fft = fftw_plan_dft_r2c_2d(_nrow, _ncol, in, inter, FFTW_MEASURE);
  rfft = fftw_plan_dft_c2r_2d(_nrow, _ncol, inter, out, FFTW_MEASURE);
  fftw_free(in);
  fftw_free(out);
}

AutoCorr2d::~AutoCorr2d() {
  fftw_destroy_plan(fft);
  fftw_destroy_plan(rfft);
  fftw_free(inter);
}

void AutoCorr2d::autocorr2(double * in_new, double * out_new) {
  fftw_execute_dft_r2c(fft, in_new, inter);
  for (int i = 0; i < ntot_complex; i++) {
    inter[i][0] = inter[i][0] * inter[i][0] + inter[i][1] * inter[i][1];
    inter[i][1] = 0.0;
  }
  fftw_execute_dft_c2r(rfft, inter, out_new);
  for (size_t i = 0; i < ntot; i++) {
    out_new[i] /= ntot;
  }
}


SpatialCorr2d::SpatialCorr2d(int ncols0, int nrows0, double Lx, double Ly) :
                             ncols(ncols0), nrows(nrows0),
                             AutoCorr2d(nrows0, ncols0) {
  ncells = ncols * nrows;
  lx = Lx / ncols;
  ly = Ly / nrows;
  dA = lx * ly;
}

void SpatialCorr2d::eval(const int *num, const double *vx, const double *vy,
                         double *c_rho, double *c_v,
                         double &rho_m, double &vx_m, double &vy_m) {
  rho_m = 0;        
  vx_m = vy_m = 0;  // velocity averaged over all particles
  double *rho = (double *)fftw_malloc(sizeof(double) * ncells);
  for (size_t i = 0; i < ncells; i++) {
    rho[i] = num[i] / dA;
    rho_m += rho[i];
  }
  autocorr2(rho, c_rho);

  double *v_tmp = (double *)fftw_malloc(sizeof(double) * ncells);
  for (size_t i = 0; i < ncells; i++) {
    v_tmp[i] = vx[i] * rho[i];
    vx_m += v_tmp[i] * dA;
  }
  double *c_vx = (double *)fftw_malloc(sizeof(double) * ncells);
  autocorr2(v_tmp, c_vx);

  for (size_t i = 0; i < ncells; i++) {
    v_tmp[i] = vy[i] * rho[i];
    vy_m += v_tmp[i] * dA;
  }
  double *c_vy = (double *)fftw_malloc(sizeof(double) * ncells);
  autocorr2(v_tmp, c_vy);

  rho_m /= ncells;
  vx_m /= (rho_m * ncells);
  vy_m /= (rho_m * ncells);
  for (size_t i = 0; i < ncells; i++) {
    c_rho[i] /= ncells;
    c_v[i] = (c_vx[i] + c_vy[i]) / (c_rho[i] * ncells);
  }

  fftw_free(rho);
  fftw_free(v_tmp);
  fftw_free(c_vx);
  fftw_free(c_vy);
}

CircleAverage::CircleAverage(int _ncols, int _nrows, double l, vector<double> &r_arr):
                             ncols(_ncols), nrows(_nrows) {
  half_nrows = nrows / 2;
  half_ncols = ncols / 2;
  RR_max = half_nrows * half_nrows;
  std::map<unsigned int, unsigned int> dict_count;
  for (int row = 0; row < half_nrows; row++) {
    unsigned int yy = row * row;
    if (row == 0) {
      for (int col = 0; col < half_ncols; col++) {
        unsigned int  RR = col * col + yy;
        if (RR < RR_max) {
          dict_count[RR] += 1;
        }
      }
    } else {
    for (int col = 0; col < ncols; col++) {
      int dx = col < half_ncols ? col : ncols - col;
      unsigned int RR = dx * dx + yy;
      if (RR < RR_max) {
        dict_count[RR] += 1;
      }
    }
    }
  }
  //r.reserve(dict_count.size());
  //count.reserve(dict_count.size());

  R_thresh = l * 10;
  RR_thresh = R_thresh * R_thresh;
  int r_ceiling = R_thresh;
  int rr_ceiling = RR_thresh;
  auto it = dict_count.cbegin();
  while (it != dict_count.cend()) {
    if (it->first < RR_thresh) {
      r_arr.push_back(sqrt(it->first));
      count.push_back(it->second);
      ++it;
    } else {
      if (rr_ceiling == RR_thresh)
        idx_thresh = r_arr.size();
      rr_ceiling += 2 * r_ceiling + 1;
      r_ceiling++;
      double sum_r = 0;
      int sum_count = 0;
      while (it->first < rr_ceiling && it != dict_count.cend()) {
        sum_r += sqrt(it->first) * it->second;
        sum_count += it->second;
        ++it;
      }
      r_arr.push_back(sum_r / sum_count);
      count.push_back(sum_count);
    }
  }
  size_r = r_arr.size();
  for (int i = 0; i < size_r; i++) {
    r_arr[i] *= l;
  }
}

void CircleAverage::eval(const double *corr2d, double *corr_r) const {
  for (int i = 0; i < size_r; i++) {
    corr_r[i] = 0;
  }
  map<unsigned int, double> corr1d;
  for (int row = 0; row < half_nrows; row++) {
    if (row == 0) {
      for (int col = 0; col < half_ncols; col++) {
        unsigned int rr = col * col;
        if (rr < RR_thresh) {
          corr1d[rr] += corr2d[col];
        } else if (rr < RR_max) {
          int idx = idx_thresh + col - R_thresh;
          corr_r[idx] += corr2d[col];
        }
      }
    } else {
      unsigned int yy = row * row;
      for (int col = 0; col < ncols; col++) {
        int x = col < half_ncols ? col : ncols - col;
        unsigned int rr = yy + x * x;
        if (rr < RR_thresh) {
          corr1d[rr] += corr2d[col + ncols * row];
        } else if (rr < RR_max) {
          int idx = idx_thresh + int(sqrt(rr)) - R_thresh;
          corr_r[idx] += corr2d[col + ncols * row];
        }
      }
    }
  }
  int j = 0;
  for (auto it = corr1d.cbegin(); it != corr1d.cend(); ++it) {
    corr_r[j] = it->second / count[j];
    j++;
  }
  for (; j < size_r; j++) {
    corr_r[j] /= count[j];
  }
}