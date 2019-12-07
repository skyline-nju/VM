#pragma once
#include "vect.h"
#include "comn.h"
#define RAND_SEQ
#define SNAP_ON

class Par_2 {
public:
  Par_2(int s, int x, int y): spin(s), pos(x, y) {}

  template <typename TRan>
  Par_2(TRan &myran, const Vec_2<int> &l);

  template <typename TRan>
  Par_2(TRan &myran, const Vec_2<int> &l, int s,
        const Vec_2<int> &origin = Vec_2<int>());

  int spin;
  Vec_2<int> pos;
};

template<typename TRan>
Par_2::Par_2(TRan & myran, const Vec_2<int> &l) {
  if (myran.doub() < 0.5) {
    spin = 1;
  } else {
    spin = -1;
  }
  pos.x = int(l.x * myran.doub());
  pos.y = int(l.y * myran.doub());
}

template <typename TRan>
Par_2::Par_2(TRan& myran, const Vec_2<int>& l, int s,
             const Vec_2<int> &origin): spin(s) {
  pos.x = int(origin.x + l.x * myran.doub());
  pos.y = int(origin.y + l.y * myran.doub());
}

class lattice_2 {
public:
  lattice_2(const Vec_2<int> &l, double beta, double eps, double rho0);

  lattice_2(const Vec_2<int> &l, double beta, double eps,
            double rho0, double h0);

  ~lattice_2();

  template <typename TRan>
  void ini_rand(TRan &myran);

  template <typename TRan>
  void ini_rand(TRan &myran, int spin0);

  template <typename TRan>
  void ini_rand(TRan &myran, int spin0, double rho_liquid, double rho_gas);

  void set_cum_pr();

  void add_particle(const Par_2 &p) const;

  void del_particle(const Par_2 &p) const;

  void hop(Par_2 &p, double rand_val) const;

  void flip(Par_2 &p, double rand_val) const;

  void flip(Par_2 &p, double rand_val, double h) const;

  double cal_m_mean() const;

  double cal_rho0() const;

  template <typename TRan>
  void one_step(TRan &myran);

  template <typename TRan>
  void one_step(TRan &myran, double h);

  void output_snap(std::ofstream &fout, int i_step);

private:
  Vec_2<int> l_;
  int n_sites_;
  unsigned short* rho_;
  short *m_;

  double beta_;
  double eps_;
  double coeff_diffusion_ = 1;
  double delta_t_;
  double cum_pr_[4]{};

  double rho0_;
  int n_par_;
  std::vector<Par_2> p_arr_;
};

template <typename TRan>
void lattice_2::ini_rand(TRan& myran) {
  for (int i = 0; i < n_par_; i++) {
    p_arr_.emplace_back(myran, l_);
    add_particle(p_arr_.back());
  }
}

template <typename TRan>
void lattice_2::ini_rand(TRan& myran, int spin0) {
  for (int i = 0; i < n_par_; i++) {
    p_arr_.emplace_back(myran, l_, spin0);
    add_particle(p_arr_.back());
  }
}

template <typename TRan>
void lattice_2::ini_rand(TRan& myran, int spin0,
                         double rho_liquid, double rho_gas) {
  int lx1 = int((rho0_ - rho_gas) / (rho_liquid - rho_gas) * l_.x);
  int n_liquid = int(lx1 * l_.y * rho_liquid);
  int n_gas = n_par_ - n_liquid;
  for (int i = 0; i < n_liquid; i++) {
    p_arr_.emplace_back(myran, Vec_2<int>(lx1, l_.y), spin0);
    add_particle(p_arr_.back());
  }
  for (int i = 0; i < n_gas; i++) {
    p_arr_.emplace_back(myran, Vec_2<int>(l_.x - lx1, l_.y), spin0,
                        Vec_2<int>(lx1, 0));
    add_particle(p_arr_.back());
  }
}

template <typename TRan>
void lattice_2::one_step(TRan& myran) {
#ifdef RAND_SEQ
  shuffle(p_arr_, myran);
  for (int i = 0; i < n_par_; i++) {
    const double rand_val = myran.doub();
    if (rand_val < cum_pr_[3]) {
      hop(p_arr_[i], rand_val);
    } else {
      flip(p_arr_[i], rand_val);
    }
  }
#else
  for (int i = 0; i < n_par_; i++) {
    const int j = int(myran.doub() * n_par_);
    const double rand_val = myran.doub();
    if (rand_val < cum_pr_[3]) {
      hop(p_arr_[j], rand_val);
    } else {
      flip(p_arr_[j], rand_val);
    }
  }
#endif
}

template <typename TRan>
void lattice_2::one_step(TRan& myran, double h) {
#ifdef RAND_SEQ
  shuffle(p_arr_, myran);
  for (int i = 0; i < n_par_; i++) {
    const double rand_val = myran.doub();
    if (rand_val < cum_pr_[3]) {
      hop(p_arr_[i], rand_val);
    } else {
      flip(p_arr_[i], rand_val, h);
    }
  }
#else
  for (int i = 0; i < n_par_; i++) {
    const int j = int(myran.doub() * n_par_);
    const double rand_val = myran.doub();
    if (rand_val < cum_pr_[3]) {
      hop(p_arr_[j], rand_val);
    } else {
      flip(p_arr_[j], rand_val, h);
    }
  }
#endif
}

void run_osc(int L, double beta, double eps, double rho0, double h0,
             int t_half, int n_period, unsigned long long seed,
             bool flag_time_ave=true);


void run_osc(int Lx, int Ly, double beta, double eps, double rho0, double h0,
             int t_half, int n_period, unsigned long long seed,
             bool flag_time_ave = true, bool rand_ini_pos = true);

void run_reverse(int L, double beta, double eps, double rho0, double h0,
                 unsigned long long seed, int n_step);