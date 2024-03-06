#include "nematic_lattice_gas.h"
#include <cmath>
#include <fstream>
#include <iomanip>

lattice_2::lattice_2(const Vec_2<int>& l, double beta, double eps, double rho0, double D)
  : l_(l), n_sites_(l.x * l.y), beta_(beta), eps_(eps), rho0_(rho0), D_(D) {
  inverse_delta_t_ = (4 * D_ + std::exp(beta_));
  rotation_rate_threshold_ = D_ * 4;

  rho_ = new unsigned short[n_sites_] {};
  Q_xx_ = new unsigned short[n_sites_] {};
  Q_yy_ = new unsigned short[n_sites_] {};
  n_par_ = int(l_.x * l_.y * rho0_);
  p_arr_.reserve(n_par_);
  std::cout << "rho0 = " << rho0_ << "\n";
  std::cout << "L = " << l_ << "\n";
  std::cout << "eps = " << eps_ << "\n";
  std::cout << "beta = " << beta_ << "\n";
  std::cout << "total particles = " << n_par_ << std::endl;
}

lattice_2::~lattice_2() {
  delete[] rho_;
  delete[] Q_xx_;
  delete[] Q_yy_;
}

void lattice_2::add_particle(const Par_2& p) const {
  const int j = p.pos.x + p.pos.y * l_.x;
  rho_[j] += 1;
  if (p.spin.y == 0) {
    Q_xx_[j] += 1;
  } else {
    Q_yy_[j] += 1;
  }
}

void lattice_2::del_particle(const Par_2& p) const {
  const int j = p.pos.x + p.pos.y * l_.x;
  rho_[j] -= 1;
  if (p.spin.y == 0) {
    Q_xx_[j] -= 1;
  } else {
    Q_yy_[j] -= 1;
  }
}

void lattice_2::hop(Par_2 & p, double rand_val) const {
  del_particle(p);
  if (p.spin.y == 0) {
    if (rand_val < D_) {
      p.pos.y += 1;
      tangle_1(p.pos.y, l_.y);
    } else if (rand_val < D_ + D_) {
      p.pos.y -= 1;
      tangle_1(p.pos.y, l_.y);
    } else if (rand_val < D_ * (3 + eps_)) {
      p.pos.x += p.spin.x;
      tangle_1(p.pos.x, l_.x);
    } else {
      p.pos.x -= p.spin.x;
      tangle_1(p.pos.x, l_.x);
    }
  } else {
    if (rand_val < D_) {
      p.pos.x += 1;
      tangle_1(p.pos.x, l_.x);
    } else if (rand_val < D_ + D_) {
      p.pos.x -= 1;
      tangle_1(p.pos.x, l_.x);
    } else if (rand_val < D_ * (3 + eps_)) {
      p.pos.y += p.spin.y;
      tangle_1(p.pos.y, l_.y);
    } else {
      p.pos.y -= p.spin.y;
      tangle_1(p.pos.y, l_.y);
    }
  }
  add_particle(p);
}

void lattice_2::flip(Par_2& p, double rand_val) const {
  int j = p.pos.x + p.pos.y * l_.x;
  double tmp = -beta_ * (Q_xx_[j] - Q_yy_[j]) / rho_[j];
  double rate = rand_val - rotation_rate_threshold_;
  if (p.spin.y == 0) {
    double w = 0.5 * std::exp(tmp);
    if (rate < w) {
      p.spin.x = 0;
      p.spin.y = 1;
      Q_xx_[j] -= 1;
      Q_yy_[j] += 1;
    } else if (rate < w + w) {
      p.spin.x = 0;
      p.spin.y = -1;
      Q_xx_[j] -= 1;
      Q_yy_[j] += 1;
    }
  } else {
    double w = 0.5 * std::exp(-tmp);
    if (rate < w) {
      p.spin.x = 1;
      p.spin.y = 0;
      Q_xx_[j] += 1;
      Q_yy_[j] -= 1;
    } else if (rate < w + w) {
      p.spin.x = -1;
      p.spin.y = 0;
      Q_xx_[j] += 1;
      Q_yy_[j] -= 1;
    }
  }
}

void lattice_2::output_snap(std::ofstream& fout) {

  unsigned short* n_arr = new unsigned short[n_sites_ * 4] {};
  // n_sites_ rows, for each row, the particle number with spins +0, -0, 0+, 0-

  for (const auto& p : p_arr_) {
    int j = p.pos.x + p.pos.y * l_.x;
    if (p.spin.y == 0) {
      if (p.spin.x == 1) {
        n_arr[4 * j]++;
      } else {
        n_arr[4 * j + 1]++;
      }
    } else {
      if (p.spin.y == 1) {
        n_arr[4 * j + 2]++;
      } else {
        n_arr[4 * j + 3]++;
      }
    }
  }
  fout.write((const char*)&n_arr[0], sizeof(unsigned short) * n_sites_ * 4);
}

void lattice_2::cal_mean_nematic_fields(double& Q_xx_mean, double& Q_yy_mean) const{
  Q_xx_mean = Q_yy_mean = 0;

  int n_tot = 0;
  for (int i = 0; i < n_sites_; i++) {
    Q_xx_mean += Q_xx_[i];
    Q_yy_mean += Q_yy_[i];
    n_tot += rho_[i];
  }
  Q_xx_mean /= n_par_;
  Q_yy_mean /= n_par_;
}


void run(int Lx, int Ly, double beta, double eps, double rho0, double D,
         unsigned long long seed, int n_step, int dn_out) {
  Vec_2<int> l(Lx, Ly);
  Ranq2 myran(seed);
  lattice_2 domain(l, beta, eps, rho0, D);
  //domain.ini_rand(myran);
  domain.ini_ordered(myran);
  char folder[100] = "data";
  char basename[100];
  snprintf(basename, 100, "L%d_%d_b%g_e%g_r%g_D%g_s%llu",
           Lx, Ly, beta, eps, rho0, D, seed);

  char order_para_file[256];
  char snap_file[256];

  int t_start = 0;
  snprintf(order_para_file, 256, "%s/%s_t%d.dat", folder, basename, t_start);
  std::ofstream fout(order_para_file);
 
  snprintf(snap_file, 256, "%s/%s_dt%d_t%d.bin", folder, basename, dn_out, t_start);
  std::ofstream fsnap(snap_file, std::ios::binary);

  double Q_xx, Q_yy;
  domain.cal_mean_nematic_fields(Q_xx, Q_yy);
  for (int i = 1; i <= n_step; i++) {
    domain.one_step(myran);
    if (i % 100 == 0) {
      domain.cal_mean_nematic_fields(Q_xx, Q_yy);
      fout << i << "\t" << Q_xx << "\t" << Q_yy << "\n";
      if (i % dn_out) {
        domain.output_snap(fsnap);
      }
    }
  }
}
