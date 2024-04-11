#pragma once
#include "ACMq4.h"

class lattice_QS_2: public lattice_2 {
public:
  lattice_QS_2(const Vec_2<int>& l, double beta,
    double eps, double rho0, double rho_thresh, double D, double eta);

  template <typename TRan>
  void one_step(TRan& myran);

private:
  double eta_;
  double rho_threshold_;
};



template <typename TRan>
void lattice_QS_2::one_step(TRan& myran) {
  shuffle(p_arr_, myran);
  for (int i = 0; i < n_par_; i++) {
    const double rand_val = myran.doub() / delta_t_;
    if (rand_val < prob_arr_[3]) {
      const int j = (p_arr_[i].pos.x + p_arr_[i].pos.y * l_.x) * 4;
      double rho_local = double(sigma_[j]) + double(sigma_[j + 1]) + double(sigma_[j + 2]) + double(sigma_[j + 3]);
      double advance_rate = D_ + eps_ * (1 + tanh(eta_ * (rho_local - rho_threshold_) / rho_threshold_));
      if (rand_val < prob_arr_[2] + advance_rate) {
        hop(p_arr_[i], rand_val);
      }
    } else {
      rotate(p_arr_[i], rand_val);
    }
  }
}


void run_QS(int Lx, int Ly, double rho0, double beta, double eps,
  double D, double rho_thresh, double eta, int n_step, int dn_out, int seed);