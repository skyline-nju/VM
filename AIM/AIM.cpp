#include "AIM.h"
#include <cmath>
#include <fstream>
#include <iomanip>

lattice_2::lattice_2(const Vec_2<int>& l, double beta, double eps, double rho0)
  : l_(l), n_sites_(l.x * l.y), beta_(beta), eps_(eps), rho0_(rho0) {
  delta_t_ = 1 / (4 * coeff_diffusion_ + std::exp(beta_));
  set_cum_pr();
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

lattice_2::lattice_2(const Vec_2<int>& l, double beta, double eps,
                     double rho0, double h0)
  : l_(l), n_sites_(l.x * l.y), beta_(beta), eps_(eps), rho0_(rho0) {
  delta_t_ = 1 / (4 * coeff_diffusion_ + std::exp(beta_ * (1 + h0)));
  set_cum_pr();
  rho_ = new unsigned short[n_sites_] {};
  m_ = new short[n_sites_] {};
  n_par_ = int(l_.x * l_.y * rho0_);
  p_arr_.reserve(n_par_);
  std::cout << "rho0 = " << rho0_ << "\n";
  std::cout << "L = " << l_ << "\n";
  std::cout << "eps = " << eps_ << "\n";
  std::cout << "beta = " << beta_ << "\n";
  std::cout << "total particles = " << n_par_ << std::endl;
  std::cout << "h0 = " << h0 << std::endl;
}

lattice_2::~lattice_2() {
  delete[] rho_;
  delete[] m_;
}

void lattice_2::set_cum_pr() {
  cum_pr_[0] = delta_t_ * coeff_diffusion_;
  cum_pr_[1] = cum_pr_[0] * 2;
  cum_pr_[2] = cum_pr_[1] + coeff_diffusion_ * (1 + eps_) * delta_t_;
  cum_pr_[3] = cum_pr_[2] + coeff_diffusion_ * (1 - eps_) * delta_t_;
  std::cout << "delta_t = " << delta_t_ << "\n";
  std::cout << "cumulative probability: " << "\n";
  std::cout << cum_pr_[0] << "; " << cum_pr_[1] << "; " << cum_pr_[2] << "; "
            << cum_pr_[3] << std::endl;
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
  if (rand_val < cum_pr_[0]) {
    p.pos.y += 1;
    if (p.pos.y >= l_.y)
      p.pos.y = 0;
  } else if (rand_val < cum_pr_[1]) {
    p.pos.y -= 1;
    if (p.pos.y < 0)
      p.pos.y += l_.y;
  } else if (rand_val < cum_pr_[2]) {
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
  const double w = std::exp(-p.spin * beta_ * double(m_[j]) / rho_[j]) * delta_t_;
  if (rand_val < cum_pr_[3] + w) {
    p.spin = -p.spin;
    m_[j] += 2 * p.spin;
  }
}

void lattice_2::flip(Par_2& p, double rand_val, double h) const {
  const int j = p.pos.x + p.pos.y * l_.x;
  const double w = std::exp(-p.spin * beta_ * (double(m_[j]) / rho_[j] + h)) * delta_t_;
  if (rand_val < cum_pr_[3] + w) {
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

void lattice_2::output_snap(std::ofstream& fout, int i_step) {
  int t = i_step;
  fout.write((const char*)&t, sizeof(int));
  fout.write((const char*)&rho_[0], sizeof(unsigned short) * n_sites_);
  fout.write((const char*)&m_[0], sizeof(short) * n_sites_);
}

void run_osc(int L, double beta, double eps, double rho0, double h0,
             int t_half, int n_period, unsigned long long seed, bool flag_time_ave) {
  Vec_2<int> l(L, L);
  int period = t_half * 2;
  int n_step = period * n_period;
  std::cout << "tot steps = " << n_step << std::endl;
  Ranq2 myran(seed);
  lattice_2 domain(l, beta, eps, rho0, h0);
  domain.ini_rand(myran, 1);
  char outfile[100];
  snprintf(outfile, 100, "%d_%g_%g_%g_%g_%d_%llu.dat",
           L, beta, eps, rho0, h0, t_half, seed);
  std::ofstream fout(outfile);

#ifdef SNAP_ON
#ifdef _MSC_VER
  snprintf(outfile, 100, R"(D:\data\DPT\AIM\snap\%d_%g_%g_%g_%g_%d_%llu.bin)",
           L, beta, eps, rho0, h0, t_half, seed);
#else
  snprintf(outfile, 100, "%d_%g_%g_%g_%g_%d_%llu.bin",
           L, beta, eps, rho0, h0, t_half, seed);
#endif
  std::ofstream fsnap(outfile, std::ios::binary);
#endif

  const auto t_start = std::chrono::system_clock::now();
  for (int j = 1; j <= n_period; j++) {
    std::cout << "j = " << j << std::endl;
    double sum_m = 0;
    for (int i = 0; i < t_half; i++) {
      domain.one_step(myran, -h0);
      int t = (j - 1) * period + i;
      if (!flag_time_ave) {
        fout << t << "\t" << domain.cal_m_mean() << "\n";
      } else {
        sum_m += domain.cal_m_mean();
      }
#ifdef SNAP_ON
      if (j == 250 || j == 500) {
        domain.output_snap(fsnap, t);
      }
#endif
    }
    for (int i = 0; i < t_half; i++) {
      domain.one_step(myran, h0);
      int t = (j - 1) * period + i + t_half;
      if (!flag_time_ave) {
        fout << t << "\t" << domain.cal_m_mean() << "\n";
      } else {
        sum_m += domain.cal_m_mean();
      }
#ifdef SNAP_ON
      if (j == 250 || j == 500) {
        domain.output_snap(fsnap, t);
      }
#endif
    }
    if (flag_time_ave) {
      fout << std::setprecision(10) << sum_m / period;
      if (j % 1000 == 0) {
        const auto t_now = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = t_now - t_start;
        const auto dt = elapsed_seconds.count();
        const auto hour = int(dt / 3600);
        const auto min = int((dt - hour * 3600) / 60);
        const int sec = dt - hour * 3600 - min * 60;
        std::cout << "period = " << j << "\tstep = " << j * period
          << "\telapsed time = " << hour << ":" << min << ":" << sec << std::endl;
        fout << std::endl;
      } else {
        fout << "\n";
      }
    }   
  }
  fout.close();
#ifdef SNAP_ON
  fsnap.close();
#endif
}

void run_osc(int Lx, int Ly, double beta, double eps, double rho0,
             double h0, int t_half, int n_period,
             unsigned long long seed, bool flag_time_ave,
             bool rand_ini_pos) {
  Vec_2<int> l(Lx, Ly);
  int period = t_half * 2;
  int n_step = period * n_period;
  std::cout << "tot steps = " << n_step << std::endl;
  Ranq2 myran(seed);
  lattice_2 domain(l, beta, eps, rho0, h0);
  if (rand_ini_pos) {
    domain.ini_rand(myran, 1);
  } else {
    domain.ini_rand(myran, 1, 11.2, 2.4);
  }
  char outfile[100];
  if (Lx != Ly) {
    snprintf(outfile, 100, "%d_%d_%g_%g_%g_%g_%d_%llu.dat",
             Lx, Ly, beta, eps, rho0, h0, t_half, seed);
  } else {
    snprintf(outfile, 100, "%d_%g_%g_%g_%g_%d_%llu.dat", Lx, beta, eps, rho0, h0, t_half, seed);
  }
  std::ofstream fout(outfile);

#ifdef SNAP_ON
#ifdef _MSC_VER
  if (Lx != Ly) {
    snprintf(outfile, 100, R"(D:\data\DPT\AIM\snap\%d_%d_%g_%g_%g_%g_%d_%llu.bin)",
      Lx, Ly, beta, eps, rho0, h0, t_half, seed);
  } else {
    snprintf(outfile, 100, R"(D:\data\DPT\AIM\snap\%d_%g_%g_%g_%g_%d_%llu.bin)",
      Lx, beta, eps, rho0, h0, t_half, seed);
  }

#else
  if (Lx != Ly) {
    snprintf(outfile, 100, "%d_%d_%g_%g_%g_%g_%d_%llu.bin",
             Lx, Ly, beta, eps, rho0, h0, t_half, seed);
  } else {
    snprintf(outfile, 100, "%d_%g_%g_%g_%g_%d_%llu.bin",
             Lx, beta, eps, rho0, h0, t_half, seed);
  }

#endif
  std::ofstream fsnap(outfile, std::ios::binary);
#endif

  const auto t_start = std::chrono::system_clock::now();
  for (int j = 1; j <= n_period; j++) {
    double sum_m = 0;
    for (int i = 0; i < t_half; i++) {
      domain.one_step(myran, -h0);
      int t = (j - 1) * period + i;
      if (!flag_time_ave) {
        fout << t << "\t" << domain.cal_m_mean() << "\n";
      } else {
        sum_m += domain.cal_m_mean();
      }
    }
#ifdef SNAP_ON
    if ((j - 1) % 1000 == 0)
      domain.output_snap(fsnap, (j - 1) * period + t_half);
#endif
    for (int i = 0; i < t_half; i++) {
      domain.one_step(myran, h0);
      int t = (j - 1) * period + i + t_half;
      if (!flag_time_ave) {
        fout << t << "\t" << domain.cal_m_mean() << "\n";
      } else {
        sum_m += domain.cal_m_mean();
      }
    }
    if (flag_time_ave) {
      fout << std::setprecision(10) << sum_m / period;
      if (j % 1000 == 0) {
        const auto t_now = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = t_now - t_start;
        const auto dt = elapsed_seconds.count();
        const auto hour = int(dt / 3600);
        const auto min = int((dt - hour * 3600) / 60);
        const int sec = dt - hour * 3600 - min * 60;
        std::cout << "period = " << j << "\tstep = " << j * period
          << "\telapsed time = " << hour << ":" << min << ":" << sec << std::endl;
        fout << std::endl;
      } else {
        fout << "\n";
      }
    }
  }
  fout.close();
#ifdef SNAP_ON
  fsnap.close();
#endif

}

void run_reverse(int L, double beta, double eps, double rho0, double h0,
                 unsigned long long seed, int n_step) {
  Vec_2<int> l(L, L);
  Ranq2 myran(seed);
  lattice_2 domain(l, beta, eps, rho0, h0);
  domain.ini_rand(myran, 1);
  char outfile[100];
  snprintf(outfile, 100, "%d_%g_%g_%g_%g_%llu.dat",
           L, beta, eps, rho0, h0, seed);
  std::ofstream fout(outfile);
  for (int i = 0; i < n_step; i++) {
    domain.one_step(myran, h0);
    fout << i << "\t" << domain.cal_m_mean() << "\n";
  }
  for (int i = n_step; i < n_step * 2; i++) {
    domain.one_step(myran, -h0);
    fout << i << "\t" << domain.cal_m_mean() << "\n";
  }
}
