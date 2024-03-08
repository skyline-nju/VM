#pragma once
#include "vect.h"
#include "comn.h"

class Par_2 {
public:
  Par_2(int s, int x, int y): spin(s), pos(x, y) {}

  template <typename TRan>
  Par_2(TRan &myran, const Vec_2<int> &l);

  template <typename TRan>
  Par_2(TRan& myran, const Vec_2<int>& l, int s0);

  Vec_2<short> pos;
  char spin;
};

template<typename TRan>
Par_2::Par_2(TRan & myran, const Vec_2<int> &l) {
  if (myran.doub() < 0.5) {
    spin = 1;
  } else {
    spin = -1;
  }
  pos.x = short(l.x * myran.doub());
  pos.y = short(l.y * myran.doub());
}

template<typename TRan>
Par_2::Par_2(TRan& myran, const Vec_2<int>& l, int s0) {
  spin = s0;
  pos.x = short(l.x * myran.doub());
  pos.y = short(l.y * myran.doub());
}

class lattice_2 {
public:
  lattice_2(const Vec_2<int> &l, double beta,
            double eps, double rho0, double D);

  ~lattice_2();

  template <typename TRan>
  void ini_rand(TRan &myran);

  template <typename TRan>
  void ini_rand(TRan &myran, int spin0);

  void add_particle(const Par_2 &p) const;

  void del_particle(const Par_2 &p) const;

  void hop(Par_2 &p, double rand_val) const;

  void flip(Par_2 &p, double rand_val) const;

  void flip(Par_2& p, double rand_val, double alpha) const;

  double cal_m_mean() const;

  double cal_rho0() const;

  template <typename TRan>
  void one_step(TRan& myran);

  template <typename TRan>
  void one_step(TRan &myran, double alpha);


  void output_snap(std::ofstream &fout);

private:
  Vec_2<int> l_;
  int n_sites_;
  unsigned short* rho_;
  short *m_;

  double beta_;
  double eps_;
  double D_ = 1;
  double delta_t_;
  double prob_arr_[4]{};

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
void lattice_2::one_step(TRan& myran) {
  shuffle(p_arr_, myran);
  for (int i = 0; i < n_par_; i++) {
    const double rand_val = myran.doub() / delta_t_;
    if (rand_val < prob_arr_[3]) {
      hop(p_arr_[i], rand_val);
    } else {
      flip(p_arr_[i], rand_val);
    }
  }
}

template <typename TRan>
void lattice_2::one_step(TRan& myran, double alpha) {
  shuffle(p_arr_, myran);
  for (int i = 0; i < n_par_; i++) {
    const double rand_val = myran.doub() / delta_t_;
    if (rand_val < prob_arr_[3]) {
      hop(p_arr_[i], rand_val);
    } else {
      flip(p_arr_[i], rand_val, alpha);
    }
  }
}

void run(int Lx, int Ly, double rho0, double beta, double eps,
         double D, int n_step, int dn_out, int seed, double alpha);

