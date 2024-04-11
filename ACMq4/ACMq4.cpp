#include "ACMq4.h"
#include <cmath>
#include <fstream>
#include <iomanip>

lattice_2::lattice_2(const Vec_2<int>& l, double beta, double eps, double rho0, double D)
  : l_(l), n_sites_(l.x * l.y), beta_(beta), eps_(eps), rho0_(rho0), D_(D) {
  prob_arr_[0] = D_;
  prob_arr_[1] = prob_arr_[0] + D_;
  prob_arr_[2] = prob_arr_[1] + 1 * (1 - eps_);
  prob_arr_[3] = prob_arr_[2] + 1 * (1 + eps_);
  //prob_arr_[2] = prob_arr_[1] + eps_;
  //prob_arr_[3] = prob_arr_[2];
  delta_t_ = 1 / (prob_arr_[3] + 2 * w0 * std::exp(beta_ * sqrt(2.) * 0.5));


  sigma_ = new unsigned char[n_sites_ * 4] {};
  n_par_ = int(l_.x * l_.y * rho0_);
  p_arr_.reserve(n_par_);
  std::cout << "rho0 = " << rho0_ << "\n";
  std::cout << "L = " << l_ << "\n";
  std::cout << "eps = " << eps_ << "\n";
  std::cout << "beta = " << beta_ << "\n";
  std::cout << "total particles = " << n_par_ << std::endl;
}

lattice_2::~lattice_2() {
  delete []sigma_;
}


void lattice_2::add_particle(const Par_2& p) const {
  const int j = (p.pos.x + p.pos.y * l_.x) * 4 + p.spin;
  sigma_[j] += 1;
}

void lattice_2::del_particle(const Par_2& p) const {
  const int j = (p.pos.x + p.pos.y * l_.x) * 4 + p.spin;
  sigma_[j] -= 1;
}

void lattice_2::hop(Par_2 & p, double rand_val) const {
  del_particle(p);
  if (p.spin == 0 || p.spin == 2) {
    if (rand_val < prob_arr_[0]) {
      p.pos.y += 1;
      if (p.pos.y >= l_.y)
        p.pos.y = 0;
    } else if (rand_val < prob_arr_[1]) {
      p.pos.y -= 1;
      if (p.pos.y < 0)
        p.pos.y += l_.y;
    } else {
      int dx = 1 - p.spin;
      if (rand_val < prob_arr_[2]) {
        p.pos.x -= dx;
        tangle_1(p.pos.x, l_.x);
      } else {
        p.pos.x += dx;
        tangle_1(p.pos.x, l_.x);
      }
    }
  } else {  // spin == 1 or 3
    if (rand_val < prob_arr_[0]) {
      p.pos.x += 1;
      if (p.pos.x >= l_.x)
        p.pos.x = 0;
    } else if (rand_val < prob_arr_[1]) {
      p.pos.x -= 1;
      if (p.pos.x < 0)
        p.pos.x += l_.x;
    } else {
      int dy = 2 - p.spin;
      if (rand_val < prob_arr_[2]) {
        p.pos.y -= dy;
        tangle_1(p.pos.y, l_.y);
      } else {
        p.pos.y += dy;
        tangle_1(p.pos.y, l_.y);
      }
    }
  }
  add_particle(p);
}

void lattice_2::rotate(Par_2& p, double rand_val) const {
  const int j = (p.pos.x + p.pos.y * l_.x) * 4;

  int rho = int(sigma_[j]) + sigma_[j + 1] + sigma_[j + 2] + sigma_[j + 3];
  int mx = int(sigma_[j]) - int(sigma_[j + 2]);
  int my = int(sigma_[j + 1]) - int(sigma_[j + 3]);

  Vec_2<int> u = get_ori(p.spin);
  char next_spin = p.spin + 1;
  if (next_spin == 4)
    next_spin = 0;
  Vec_2<int> next_u = get_ori(next_spin);
  Vec_2<int> du = next_u - u;
  
  const double w_next = w0 * std::exp(0.5 * beta_ / rho * (mx * du.x + my * du.y));
  if (rand_val < prob_arr_[3] + w_next) {
    sigma_[j + p.spin] -= 1;
    p.spin = next_spin;
    sigma_[j + p.spin] += 1;
  } else {
    char prev_spin = p.spin - 1;
    if (prev_spin == -1) {
      prev_spin = 3;
    }
    Vec_2<int> prev_u = get_ori(prev_spin);
    du = prev_u - u;
    const double w_prev = w0 * std::exp(0.5 * beta_ / rho * (mx * du.x + my * du.y));
    if (rand_val < prob_arr_[3] + w_next + w_prev) {
      sigma_[j + p.spin] -= 1;
      p.spin = prev_spin;
      sigma_[j + p.spin] += 1;
    }
  }
}


void lattice_2::cal_m_mean(double &mx, double &my) const {
  mx = 0;
  my = 0;
  for (int i = 0; i < n_sites_; i++) {
    mx += (int(sigma_[i * 4]) - sigma_[i * 4 + 2]);
    my += (int(sigma_[i * 4 + 1]) - sigma_[i * 4 + 3]);
  }
  mx /= n_par_;
  my /= n_par_;
}

double lattice_2::cal_rho0() const {
  double rho_sum = 0;
  for (int i = 0; i < n_sites_; i++) {
    rho_sum += (int(sigma_[i * 4]) + sigma_[i * 4 + 1] + sigma_[i * 4 + 2] + sigma_[i * 4 + 3]);
  }
  return rho_sum / n_sites_;
}

void lattice_2::output_snap(std::ofstream& fout) {
  fout.write((const char*)&sigma_[0], sizeof(unsigned char) * n_sites_ * 4);
}

void run(int Lx, int Ly, double rho0, double beta, double eps,
         double D, int n_step, int dn_out, int seed) {
  Vec_2<int> l(Lx, Ly);
  Ranq2 myran(seed);
  lattice_2 domain(l, beta, eps, rho0, D);
  domain.ini_rand(myran, 0);
  //domain.ini_rand(myran);
  char folder[100] = "data";
  char basename[100];
  snprintf(basename, 100, "L%d_%d_b%g_r%g_e%g_D%g_s%d",
    Lx, Ly, beta, rho0, eps, D, seed);

  char order_para_file[256];
  char snap_file[256];

  int t_start = 0;
  snprintf(order_para_file, 256, "%s/%s_t%d.dat", folder, basename, t_start);
  std::ofstream fout(order_para_file);

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



