#pragma once
#include "vect.h"
#include "comn.h"
#define RAND_SEQ
#define SNAP_ON

class Par_2 {
public:
  Par_2(int s_x, int s_y, int x, int y): spin(s_x, s_y), pos(x, y) {}

  template <typename TRan>
  Par_2(TRan &myran, const Vec_2<int> &l);

  template <typename TRan>
  Par_2(TRan &myran, const Vec_2<int> &l, const Vec_2<short> &s,
        const Vec_2<int> &origin = Vec_2<int>());

  Vec_2<short> spin;
  Vec_2<int> pos;
};

template<typename TRan>
Par_2::Par_2(TRan & myran, const Vec_2<int> &l) {
  double p = myran.doub();
  if (p < 0.25) {
    spin = Vec_2<short>(1, 0);
  } else if (p < 0.5) {
    spin = Vec_2<short>(-1, 0);
  } else if (p < 0.75) {
    spin = Vec_2<short>(0, 1);
  } else {
    spin = Vec_2<short>(0, -1);
  }
  pos.x = int(l.x * myran.doub());
  pos.y = int(l.y * myran.doub());
}

template <typename TRan>
Par_2::Par_2(TRan& myran, const Vec_2<int>& l, const Vec_2<short>& s,
             const Vec_2<int> &origin): spin(s) {
  pos.x = int(origin.x + l.x * myran.doub());
  pos.y = int(origin.y + l.y * myran.doub());
}

class lattice_2 {
public:
  lattice_2(const Vec_2<int> &l, double beta, double eps, double rho0, double D);

  ~lattice_2();

  template <typename TRan>
  void ini_rand(TRan &myran);

  template <typename TRan>
  void ini_rand(TRan &myran, Vec_2<short>& spin0);

  template <typename TRan>
  void ini_rand(TRan &myran, Vec_2<short>& spin0, double rho_liquid, double rho_gas);

  template <typename TRan>
  void ini_ordered(TRan& myran);


  void add_particle(const Par_2 &p) const;

  void del_particle(const Par_2 &p) const;

  void hop(Par_2 &p, double rand_val) const;

  void flip(Par_2 &p, double rand_val) const;

  template <typename TRan>
  void hop_or_flip(Par_2& p, TRan& myran) const;

  template <typename TRan>
  void one_step(TRan &myran);

  void cal_mean_nematic_fields(double& Q_xx_mean, double& Q_yy_mean) const;


  void output_snap(std::ofstream &fout);

private:
  Vec_2<int> l_;
  int n_sites_;
  unsigned short* rho_;
  unsigned short* Q_xx_;
  unsigned short* Q_yy_;

  double beta_;
  double eps_;
  double D_;
  double inverse_delta_t_;
  double rotation_rate_threshold_;

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
  std::cout << "created " << n_par_ << " particles with random position and orientation" << std::endl;
}

template <typename TRan>
void lattice_2::ini_rand(TRan& myran, Vec_2<short>& spin0) {
  for (int i = 0; i < n_par_; i++) {
    p_arr_.emplace_back(myran, l_, spin0);
    add_particle(p_arr_.back());
  }
}

template <typename TRan>
void lattice_2::ini_rand(TRan& myran, Vec_2<short>& spin0,
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

template<typename TRan>
void lattice_2::ini_ordered(TRan& myran) {
  for (int i = 0; i < n_par_; i++) {
    p_arr_.emplace_back(myran, l_);
    if (myran.doub() > 0.5) {
      p_arr_[i].spin.x = 1;
    } else {
      p_arr_[i].spin.x = -1;
    }
    p_arr_[i].spin.y = 0;
    add_particle(p_arr_[i]);
  }
  std::cout << "created " << n_par_ << " particles with random position but ordered orientation" << std::endl;
}

template<typename TRan>
void lattice_2::hop_or_flip(Par_2& p, TRan& myran) const {
  double rand_val = myran.doub() * inverse_delta_t_;
  if (rand_val < rotation_rate_threshold_) {
    hop(p, rand_val);
  } else {
    flip(p, rand_val);
  }
}

template <typename TRan>
void lattice_2::one_step(TRan& myran) {
  shuffle(p_arr_, myran);
  for (int i = 0; i < n_par_; i++) {
    hop_or_flip(p_arr_[i], myran);
  }
}


void run(int Lx, int Ly, double beta, double eps, double rho0, double D,
  unsigned long long seed, int n_step, int dn_out);