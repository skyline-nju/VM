#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include "particle.h"
#include "rand.h"
#include "comn.h"
#include "grid.h"

// Base class for Vicsek model
class VM {
public:
  VM(double Lx, double Ly, int N, double eta, double v0=0.5)
  : Lx_(Lx), Ly_(Ly), half_Lx_(Lx * 0.5), half_Ly_(Ly * 0.5), N_(N), eta_(eta), v0_(v0) {}

  ~VM(){};
 
  // get data
  int get_n_birds() const { return N_; }

protected:
  double v0_;
  double Lx_;
  double Ly_;
  double half_Lx_;
  double half_Ly_;
  int N_;
  double eta_;
};

/*************************************************************************************
 *      Class for one-species flocks with metric interaction
 * 
 * 
 ************************************************************************************/
template <class TNode>
class VM_1: public VM {
public:
  VM_1(double Lx, double Ly, int N, double eta, double v0=0.5);

  template <typename T>
  void input_data(const T *x, const T *y, const T *vx, const T *vy);

  template <typename T>
  void get_pos_arr(T* pos) const;

  template <typename TRan>
  void ini_rand(TRan &myran);

  template <typename TRan>
  void ini_rand(TRan &myran, double theta0);

  template <typename TSnap>
  void ini_from_snap(TSnap& snap_reader);

  template <typename TSnap>
  void dump(int i_step, TSnap& snap_writer);

  void align();

  template<typename TRan>
  void stream(TRan& myran);

  void get_x(int i, double &x, double &y) const;
  void get_v(int i, double &vx, double &vy) const { vx = p_arr_[i].vx; vy = p_arr_[i].vy; }
  void get_theta(int i, double & Theta) const { Theta = p_arr_[i].get_theta(); }
  template <typename T>
  void get_order_para(T &phi, T &theta) const;

protected:
  // linked node celllist
  Grid<TNode> cl_;

  // array of particles
  std::vector<TNode> p_arr_;

};

template <class TNode>
VM_1<TNode>::VM_1(double Lx, double Ly, int N, double eta, double v0)
  : VM(Lx, Ly, N, eta, v0), cl_(Lx, Ly, 1.) {
  p_arr_.reserve(N);
}

template<class TNode>
template<typename T>
void VM_1<TNode>::input_data(const T *x, const T *y, const T *vx, const T *vy) {
  for (int i = 0; i < N_; i++) {
    p_arr_.emplace_back(x[i], y[i], vx[i], vy[i]);
  }
  cl_.link_nodes(p_arr_);
}

template <class TNode>
template <typename T>
void VM_1<TNode>::get_pos_arr(T *pos) const {
  for (int i = 0; i < N_; i++) {
    size_t j = i * 3;
    double x, y, theta;
    get_x(i, x, y);
    get_theta(i, theta);
    pos[j] = x - half_Lx_;
    pos[j+1] = y - half_Ly_;
    pos[j+2] = theta;
  }
}

template <class TNode>
template <typename TRan>
void VM_1<TNode>::ini_rand(TRan &myran) {
  double *x = new double[N_];
  double *y = new double[N_];
  double *vx = new double[N_];
  double *vy = new double[N_];
  for (int i = 0; i < N_; i++) {
    x[i] = myran.doub() * Lx_;
    y[i] = myran.doub() * Ly_;
    double theta = PI * myran.doub() * 2.;
    vx[i] = cos(theta);
    vy[i] = sin(theta);
  }
  input_data(x, y, vx, vy);
  delete []x;
  delete []y;
  delete []vx;
  delete []vy;
}

template <class TNode>
template <typename TRan>
void VM_1<TNode>::ini_rand(TRan &myran, double theta0) {
  double *x = new double[N_];
  double *y = new double[N_];
  double *vx = new double[N_];
  double *vy = new double[N_];

  double vx0 = cos(theta0);
  double vy0 = sin(theta0);
  for (int i = 0; i < N_; i++) {
    x[i] = myran.doub() * Lx_;
    y[i] = myran.doub() * Ly_;
    vx[i] = vx0;
    vy[i] = vy0;
  }
  input_data(x, y, vx, vy);
  delete []x;
  delete []y;
  delete []vx;
  delete []vy;
}

template <class TNode>
template <typename TSnap>
void VM_1<TNode>::ini_from_snap(TSnap &snap_reader) {
  double *x = new double[N_];
  double *y = new double[N_];
  double *vx = new double[N_];
  double *vy = new double[N_];
  snap_reader.read_last_frame(x, y, vx, vy);
  input_data(x, y, vx, vy);
  delete []x;
  delete []y;
  delete []vx;
  delete []vy;
}

template <class TNode>
template <typename TSnap>
void VM_1<TNode>::dump(int i_step, TSnap &snap_writer) {
  if (snap_writer.need_export(i_step)) {
    float* pos = new float[N_*3];
    get_pos_arr(pos);
    snap_writer.dump(i_step, N_, pos);
    delete [] pos;
  }
}

template <class TNode>
void VM_1<TNode>::align() {
  cl_.all_pairs();
}

template <class TNode>
template <typename TRan>
void VM_1<TNode>::stream(TRan &myran) {
  for (int i = 0; i < N_; i++) {
    p_arr_[i].move(eta_, myran, Lx_, Ly_, v0_);
  }
  cl_.refresh(p_arr_);
}

template<class TNode>
void VM_1<TNode>::get_x(int i, double & x, double & y) const {
  x = p_arr_[i].x;
  y = p_arr_[i].y;
}

template <class TNode>
template <typename T>
void VM_1<TNode>::get_order_para(T &phi, T &theta) const {
  double vx = 0;
  double vy = 0;
  for (int i = 0; i < N_; i++) {
    vx += p_arr_[i].vx;
    vy += p_arr_[i].vy;
  }
  vx /= N_;
  vy /= N_;
  phi = sqrt(vx * vx + vy * vy);
  theta = atan2(vy, vx); 
}


template <class TNode>
class VM_2: public VM_1<TNode> {
public:
  VM_2(double Lx, double Ly, int N,  int n_dis, double eta, double v0=0.5):
    VM_1<TNode>(Lx, Ly, N, eta, v0), n_dis_(n_dis) {}

  template <typename TSnap>
  void dump(int i_step, TSnap& snap_writer);

  void set_type();

protected:
  int n_dis_;
};


template <class TNode>
template <typename TSnap>
void VM_2<TNode>::dump(int i_step, TSnap &snap_writer) {
  if (snap_writer.need_export(i_step)) {
    float* pos = new float[VM::N_*3];
    uint32_t* type_id = new uint32_t[VM::N_]{};
    for (int i = 0; i < n_dis_; i++) {
      type_id[i] = 1;
    }
    VM_1<TNode>::get_pos_arr(pos);
    snap_writer.dump(i_step, VM::N_, pos, type_id);
    delete [] pos;
    delete [] type_id;
  }
}

template <class TNode>
void VM_2<TNode>::set_type() {
  for (int i=0; i< VM::N_; i++) {
    if (i < n_dis_) {
      VM_1<TNode>::p_arr_[i].set_type(1);
    } else {
      VM_1<TNode>::p_arr_[i].set_type(0);
    }

  }
}