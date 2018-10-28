#ifndef EXPORTER_H
#define EXPORTER_H
#include "comn.h"
#include "cmdline.h"
#include "spatial_corr.h"
#include <iomanip>

/*************************************************************************//**
 * \brief Log exporter
 *
 ****************************************************************************/
class LogExporter :public BaseLogExporter {
public:
  explicit LogExporter(const cmdline::parser &cmd);
};

/************************************************************************//**
 * \brief Exporter for 2d spatial correlation functions
 *
 ***************************************************************************/
class OrderParaExpoter: public BaseExporter {
public:
  explicit OrderParaExpoter(const std::string &filename,
                            int end, int sep, int start = 0);
  template <typename TPar>
  void write(int i_step, const TPar *p_arr, int n_par);
private:
  std::ofstream fout_;
};

template <typename TPar>
void OrderParaExpoter::write(int i_step, const TPar* p_arr, int n_par) {
  double svx = 0;
  double svy = 0;
  for (int i = 0; i < n_par; i++) {
    svx += p_arr[i].vx;
    svy += p_arr[i].vy;
  }
  double phi = std::sqrt(svx * svx + svy * svy) / n_par;
  double theta = std::atan2(svy, svx);
  fout_ << i_step << "\t" << std::setprecision(8) << phi << "\t"
    << theta << "\n";
}

/*************************************************************************//**
 * \brief Exporter for coarse-grained snapshots
 * 
 ****************************************************************************/
class CoarseGrainSnapExporter: public BaseExporter {
public:
  CoarseGrainSnapExporter(const std::string &filename, int ncols, int nrows,
                         int end, int sep, int start = 0);

  template <typename TPar>
  void write_frame(int i_step, const TPar *p_arr, int n_par);

private:
  void write_frame(int i_step, const short *num,
                   const float *vx, const float *vy);
  int ncid_;
  int time_id_;
  int num_id_;
  int vx_id_;
  int vy_id_;

  size_t frame_len_;
  size_t time_idx_[1];

  int ncols_;
  int nrows_;
  double lx_;
  double ly_;
};

template <typename TPar>
void CoarseGrainSnapExporter::write_frame(int i_step, const TPar* p_arr,
                                         int n_par) {
  size_t ncells = ncols_ * nrows_;
  short *num = new short[ncells]{};
  float *vx = new float[ncells]{};
  float *vy = new float[ncells]{};

  for (int i = 0; i < n_par; i++) {
    int col = int(p_arr[i].x / lx_);
    int row = int(p_arr[i].y / ly_);
    int i_cell = col + row * ncols_;
    num[i_cell] += 1;
    vx[i_cell] += p_arr[i].vx;
    vy[i_cell] += p_arr[i].vy;
  }
  write_frame(i_step, num, vx, vy);
  delete[] num;
  delete[] vx;
  delete[] vy;
}
/************************************************************************//**
 * \brief Exporter for 2d spatial correlation functions
 * 
 ***************************************************************************/

class SpatialCorr2dExporter : public SpatialCorr2d, public BaseExporter {
public:
  SpatialCorr2dExporter(const std::string &filename, int ncols, int nrows,
                        int end, int sep, int start = 0);

  template <typename TPar>
  void write_frame(int i_step, const TPar *p_arr, int n_par);

private:
  void write_frame(int i_step, const float *c_rho, const float *c_v,
                   const double *v_mean);
  int ncid_;
  int time_id_;
  int corr_rho_id_;
  int corr_v_id_;
  int v_mean_id_;

  size_t frame_len_;
  size_t time_idx_[1];
};


template<typename TPar>
void SpatialCorr2dExporter::write_frame(int i_step, const TPar * p_arr, int n_par) {
  double v_mean[2] = {0., 0.};
  float *c_rho = new float[nrows_ / 2 * ncols_];
  float *c_v = new float[nrows_ / 2 * ncols_];
  cal_corr(p_arr, n_par, c_rho, c_v, v_mean[0], v_mean[1]);
  v_mean[0] /= n_par;
  v_mean[1] /= n_par;
  write_frame(i_step, c_rho, c_v, v_mean);
  delete[] c_rho;
  delete[] c_v;
}

void set_output(const cmdline::parser &cmd,
                LogExporter **log_ex,
                OrderParaExpoter **order_ex,
                CoarseGrainSnapExporter **cg_ex,
                SpatialCorr2dExporter **corr_ex);

template <typename TPar>
void output(int i_step, const TPar *p_arr, int n_par,
            LogExporter * log_ex, OrderParaExpoter *order_ex,
            CoarseGrainSnapExporter *cg_ex, SpatialCorr2dExporter *corr_ex) {
  if (log_ex) {
    log_ex->record(i_step);
  }
  if (order_ex && order_ex->need_export(i_step)) {
    order_ex->write(i_step, p_arr, n_par);
  }
  if(cg_ex && cg_ex->need_export(i_step)) {
    cg_ex->write_frame(i_step, p_arr, n_par);
  }
  if(corr_ex && corr_ex->need_export(i_step)) {
    corr_ex->write_frame(i_step, p_arr, n_par);
  }
}
#endif
