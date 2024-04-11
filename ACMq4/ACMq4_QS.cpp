#include "ACMq4_QS.h"
#include <cmath>
#include <fstream>
#include <iomanip>

lattice_QS_2::lattice_QS_2(const Vec_2<int>& l, double beta, double eps, double rho0, double rho_thresh, double D, double eta)
  : lattice_2(l, beta, eps, rho0, D), rho_threshold_(rho_thresh), eta_(eta) {
  prob_arr_[0] = D_;
  prob_arr_[1] = prob_arr_[0] + D_;
  prob_arr_[2] = prob_arr_[1] + D_;
  prob_arr_[3] = prob_arr_[2] + D_ + 2 * eps_;
  delta_t_ = 1 / (prob_arr_[3] + 2 * w0 * std::exp(beta_ * sqrt(2.) * 0.5));
}

void run_QS(int Lx, int Ly, double rho0, double beta, double eps, double D, double rho_thresh, double eta, int n_step, int dn_out, int seed) {
  Vec_2<int> l(Lx, Ly);
  Ranq2 myran(seed);
  lattice_QS_2 domain(l, beta, eps, rho0, rho_thresh, D, eta);
  domain.ini_rand(myran, 0);
  //domain.ini_rand(myran);
  char folder[100] = "data";
  char basename[100];
  snprintf(basename, 100, "L%d_%d_b%g_e%g_r%g_%g_v%g_D%g_s%d",
    Lx, Ly, beta, eta, rho0, rho_thresh, eps, D, seed);

  char order_para_file[256];
  char snap_file[256];

  int t_start = 0;
  //snprintf(order_para_file, 256, "%s/%s_t%d.dat", folder, basename, t_start);
  //std::ofstream fout(order_para_file);

  snprintf(snap_file, 256, "%s/%s_dt%d_t%d.bin", folder, basename, dn_out, t_start);
  std::ofstream fsnap(snap_file, std::ios::binary);


  for (int i = 1; i <= n_step; i++) {
    domain.one_step(myran);
    if (i % 1000 == 0) {
      double mx, my;
      domain.cal_m_mean(mx, my);
      double rho0 = domain.cal_rho0();
      std::cout << i << "\t" << rho0 << "\t" << mx << "\t" << my << std::endl;
      //fout << i << "\t" << domain.cal_m_mean() << "\n";
      if (i % dn_out == 0) {
        domain.output_snap(fsnap);
      }
    }
  }
}
