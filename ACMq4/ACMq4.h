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
  spin = char(4 * myran.doub());
  pos.x = short(l.x * myran.doub());
  pos.y = short(l.y * myran.doub());
}

template<typename TRan>
Par_2::Par_2(TRan& myran, const Vec_2<int>& l, int s0) {
  spin = s0;
  pos.x = short(l.x * myran.doub());
  pos.y = short(l.y * myran.doub());
}



template <typename TInt>
Vec_2<int> get_ori(TInt spin){
  if (spin == 0) {
    return Vec_2<int>(1, 0);
  } else if (spin == 1) {
    return Vec_2<int>(0, 1);
  } else if (spin == 2) {
    return Vec_2<int>(-1, 0);
  } else {
    return Vec_2<int>(0, -1);
  }
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

  void rotate(Par_2 &p, double rand_val) const;


  void cal_m_mean(double &mx, double &my) const;

  double cal_rho0() const;

  template <typename TRan>
  void one_step(TRan& myran);



  void output_snap(std::ofstream &fout);

protected:
  Vec_2<int> l_;
  int n_sites_;
  //unsigned short* rho_;
  //short *m_;
  unsigned char* sigma_;


  double beta_;
  double eps_;
  double D_;
  double delta_t_;
  double prob_arr_[4]{};

  double rho0_;
  int n_par_;
  std::vector<Par_2> p_arr_;
  double w0 = 4 * 4 / ( 4 * PI * PI); // q^2 / 4 PI^2
};

template <typename TRan>
void lattice_2::ini_rand(TRan& myran) {
  for (int i = 0; i < n_par_; i++) {
    p_arr_.emplace_back(myran, l_);
    add_particle(p_arr_.back());
  }
  double rho0 = cal_rho0();
  std::cout << "create " << n_par_ << " particles with density " << rho0 << std::endl;
}

template <typename TRan>
void lattice_2::ini_rand(TRan& myran, int spin0) {
  for (int i = 0; i < n_par_; i++) {
    p_arr_.emplace_back(myran, l_, spin0);
    add_particle(p_arr_.back());
  }
  double rho0 = cal_rho0();
  std::cout << "create " << n_par_ << " particles with density " << rho0 << std::endl;
}


template <typename TRan>
void lattice_2::one_step(TRan& myran) {
  shuffle(p_arr_, myran);
  for (int i = 0; i < n_par_; i++) {
    const double rand_val = myran.doub() / delta_t_;
    if (rand_val < prob_arr_[3]) {
      hop(p_arr_[i], rand_val);
    } else {
      rotate(p_arr_[i], rand_val);
    }
  }
}

void run(int Lx, int Ly, double rho0, double beta, double eps,
         double D, int n_step, int dn_out, int seed);

