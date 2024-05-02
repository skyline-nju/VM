#include "AIM.h"
#include <cmath>
#include <fstream>
#include <iomanip>

lattice_2::lattice_2(const Vec_2<int>& l, double beta, double eps, double rho0, double D)
  : l_(l), n_sites_(l.x * l.y), beta_(beta), eps_(eps), rho0_(rho0) {
  prob_arr_[0] = D_;
  prob_arr_[1] = prob_arr_[0] + D_;
  prob_arr_[2] = prob_arr_[1] + D_ * (1 + eps_);
  prob_arr_[3] = prob_arr_[2] + D_ * (1 - eps_);
  //prob_arr_[2] = prob_arr_[1] + eps_;
  //prob_arr_[3] = prob_arr_[2];
  delta_t_ = 1 / (prob_arr_[3] + std::exp(beta_));
  
  rho_ = new unsigned short[n_sites_] {};
  m_ = new short[n_sites_] {};
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
  delete[] m_;
}


void lattice_2::add_particle(const Par_2& p) const {
  const int j = p.pos.x + p.pos.y * l_.x;
  rho_[j] += 1;
  m_[j] += p.spin;
}

void lattice_2::del_particle(const Par_2& p) const {
  const int j = p.pos.x + p.pos.y * l_.x;
  rho_[j] -= 1;
  m_[j] -= p.spin;
}

void lattice_2::hop(Par_2 & p, double rand_val) const {
  del_particle(p);
  if (rand_val < prob_arr_[0]) {
    p.pos.y += 1;
    if (p.pos.y >= l_.y)
      p.pos.y = 0;
  } else if (rand_val < prob_arr_[1]) {
    p.pos.y -= 1;
    if (p.pos.y < 0)
      p.pos.y += l_.y;
  } else if (rand_val < prob_arr_[2]) {
    p.pos.x += p.spin;
    tangle_1(p.pos.x, l_.x);
  } else {
    p.pos.x -= p.spin;
    tangle_1(p.pos.x, l_.x);
  }
  add_particle(p);
}

void lattice_2::flip(Par_2& p, double rand_val) const {
  const int j = p.pos.x + p.pos.y * l_.x;
  const double w = std::exp(-p.spin * beta_ * double(m_[j]) / rho_[j]);
  if (rand_val < prob_arr_[3] + w) {
    p.spin = -p.spin;
    m_[j] += 2 * p.spin;
  }
}

void lattice_2::flip(Par_2& p, double rand_val, double alpha) const {
  int col = p.pos.x;
  int row_Lx = p.pos.y * l_.x;

  int col_next = col + p.spin;
  tangle_1(col_next, l_.x);

  int j = col + row_Lx;
  int j_next = col_next + row_Lx;

  double m_mean = ((1 - alpha) * m_[j] + alpha * m_[j_next]) / ((1 - alpha) * rho_[j] + alpha * rho_[j_next]);
  double w = std::exp(-p.spin * beta_ * m_mean);
  if (rand_val < prob_arr_[3] + w) {
    p.spin = -p.spin;
    m_[j] += 2 * p.spin;
  }
}

double lattice_2::cal_m_mean() const {
  double m_sum = 0;
  for (int i = 0; i < n_sites_; i++) {
    m_sum += m_[i];
  }
  return m_sum / n_par_;
}

double lattice_2::cal_rho0() const {
  double rho_sum = 0;
  for (int i = 0; i < n_sites_; i++) {
    rho_sum += rho_[i];
  }
  return rho_sum / n_sites_;
}

void lattice_2::output_snap(std::ofstream& fout) {
  fout.write((const char*)&rho_[0], sizeof(unsigned short) * n_sites_);
  fout.write((const char*)&m_[0], sizeof(short) * n_sites_);
}

void run(int Lx, int Ly, double rho0, double beta, double eps,
         double D, int n_step, int dn_out, int seed, double alpha,
         const std::string& ini_condi) {
  Vec_2<int> l(Lx, Ly);
  Ranq2 myran(seed);
  lattice_2 domain(l, beta, eps, rho0, D);
  char folder[100] = "data";
  char basename[100];
  snprintf(basename, 100, "L%d_%d_b%g_r%g_a%g_e%g_D%g_s%d",
    Lx, Ly, beta, rho0, alpha, eps, D, seed);

  char order_para_file[256];
  char snap_file[256];

  int t_start = 0;
  snprintf(order_para_file, 256, "%s/%s_t%d.dat", folder, basename, t_start);
  snprintf(snap_file, 256, "%s/%s_dt%d_t%d.bin", folder, basename, dn_out, t_start);

  domain.ini_particles(myran, ini_condi, snap_file);
  //std::ofstream fout(order_para_file);
  std::ofstream fsnap;

  if (ini_condi == "resume") {
    fsnap.open(snap_file, std::ios::binary | std::ios::app);
  } else {
    fsnap.open(snap_file, std::ios::binary);
  }


  for (int i = 1; i <= n_step; i++) {
    domain.one_step(myran, alpha);
    if (i % 1000 == 0) {
      std::cout << i << "\t" << domain.cal_m_mean() << "\n";
      if (i % dn_out == 0) {
        domain.output_snap(fsnap);
      }
    }
  }
}



