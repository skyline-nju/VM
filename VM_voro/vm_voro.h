#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include "particle.h"
#include "rand.h"
#include "comn.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_2.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_traits_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Vector_2.h>

// Base class for Vicsek model
class VM {
public:
  VM(double Lx, double Ly, int N, double eta, double v0=0.5)
  : Lx_(Lx), Ly_(Ly), half_Lx_(Lx * 0.5), half_Ly_(Ly * 0.5), N_(N), eta_(eta), v0_(v0) {}

  ~VM(){};
 
  // get data
  int get_num_birds() const { return N_; }

protected:
  double v0_;
  double Lx_;
  double Ly_;
  double half_Lx_;
  double half_Ly_;
  int N_;
  double eta_;
};

template <class BaseV>
class VM_voro: public VM {
public:
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef CGAL::Vector_2<K> Vec_2;
  typedef CGAL::Periodic_2_Delaunay_triangulation_traits_2<K> GT;
  typedef CGAL::Periodic_2_triangulation_vertex_base_2<GT> Vb;
  typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned int, GT, Vb> VbInfo;
  typedef CGAL::Periodic_2_triangulation_face_base_2<GT> Fb;
  typedef CGAL::Triangulation_data_structure_2<VbInfo, Fb> Tds;
  typedef CGAL::Periodic_2_Delaunay_triangulation_2<GT, Tds> PDT;
  typedef PDT::Point Point;
  typedef Tds::Vertex_handle Vertex_handle;

  VM_voro(double Lx, double Ly, int N, double eta, double v0=0.5);

  template <typename T>
  void input_data(const T *x, const T *y, const T *vx, const T *vy);

  template <typename TRan>
  void ini_rand(TRan myran);

  template <typename TRan>
  void ini_rand(TRan myran, double theta0);

  void align();

  // the first n_dis particles are dissenters
  void align(int n_dis);

  template<typename TRan>
  void stream(TRan& myran);

  void get_x(int i, double &x, double &y) const;
  void get_v(int i, double &vx, double &vy) const { v_arr_[i].get_v(vx, vy); }
  void get_theta(int i, double & Theta) const { v_arr_[i].get_theta(Theta); }
  template <typename T>
  void get_order_para(T &phi, T &theta) const;

protected:
  // array of velocity
  std::vector<BaseV> v_arr_;

  // array of vertex handle
  std::vector<Vertex_handle> vh_arr_;

  PDT *DT_;
};

template <class BaseV>
VM_voro<BaseV>::VM_voro(double Lx, double Ly, int N, double eta, double v0)
  : VM(Lx, Ly, N, eta, v0) {
  // initialize triangulation
  PDT::Iso_rectangle domain(0, 0, Lx_, Ly_);
  DT_ = new PDT(domain);
}

template<class BaseV>
template<typename T>
void VM_voro<BaseV>::input_data(const T *x, const T *y, const T *vx, const T *vy) {
  v_arr_.reserve(N_);
  vh_arr_.reserve(N_);
  for (int i = 0; i < N_; i++) {
    v_arr_.emplace_back(vx[i], vy[i]);
    auto vh = DT_->insert(Point(x[i], y[i]));
    vh->info() = i;
    vh_arr_.push_back(vh);
  }
}

template <class BaseV>
template <typename TRan>
void VM_voro<BaseV>::ini_rand(TRan myran) {
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

template <class BaseV>
template <typename TRan>
void VM_voro<BaseV>::ini_rand(TRan myran, double theta0) {
  double *x = new double[N_];
  double *y = new double[N_];
  double *vx = new double[N_];
  double *vy = new double[N_];
  for (int i = 0; i < N_; i++) {
    x[i] = myran.doub() * Lx_;
    y[i] = myran.doub() * Ly_;
    vx[i] = cos(theta0);
    vy[i] = sin(theta0);
  }
  input_data(x, y, vx, vy);
  delete []x;
  delete []y;
  delete []vx;
  delete []vy;
}

template <class BaseV>
void VM_voro<BaseV>::align() {
  for (auto fit = DT_->periodic_triangles_begin(PDT::UNIQUE);
    fit != DT_->periodic_triangles_end(PDT::UNIQUE); ++fit) {
    unsigned int idx0 = fit.get_face()->vertex(0)->info();
    unsigned int idx1 = fit.get_face()->vertex(1)->info();
    unsigned int idx2 = fit.get_face()->vertex(2)->info();
    if (idx0 < idx1) {
      v_arr_[idx0].collide(&v_arr_[idx1]);
    }
    if (idx1 < idx2) {
      v_arr_[idx1].collide(&v_arr_[idx2]);
    }
    if (idx2 < idx0) {
      v_arr_[idx2].collide(&v_arr_[idx0]);
    }
  }
}

template <class BaseV>
void VM_voro<BaseV>::align(int n_dis) {
  for (auto fit = DT_->periodic_triangles_begin(PDT::UNIQUE);
    fit != DT_->periodic_triangles_end(PDT::UNIQUE); ++fit) {
    unsigned int idx0 = fit.get_face()->vertex(0)->info();
    unsigned int idx1 = fit.get_face()->vertex(1)->info();
    unsigned int idx2 = fit.get_face()->vertex(2)->info();

    bool is_aligner_0 = idx0 >= n_dis;
    bool is_aligner_1 = idx1 >= n_dis;
    bool is_aligner_2 = idx2 >= n_dis;
    if (idx0 < idx1) {
      v_arr_[idx0].collide(&v_arr_[idx1], is_aligner_0, is_aligner_1);
    }
    if (idx1 < idx2) {
      v_arr_[idx1].collide(&v_arr_[idx2], is_aligner_1, is_aligner_2);
    }
    if (idx2 < idx0) {
      v_arr_[idx2].collide(&v_arr_[idx0], is_aligner_2, is_aligner_0);
    }
  }
}

template <class BaseV>
template <typename TRan>
void VM_voro<BaseV>::stream(TRan &myran) {
  for (int i = 0; i < N_; i++) {
    v_arr_[i].update_v(eta_, myran);
    double x, y;
    get_x(i, x, y);
    v_arr_[i].update_x(x, y, v0_, Lx_, Ly_);
    Point p_new(x, y);
    auto vh_old = vh_arr_[i];
    vh_arr_[i] = DT_->insert(p_new, vh_old->face());
    vh_arr_[i]->info() = vh_old->info();
    DT_->remove(vh_old);
  }
}

template<class BaseV>
void VM_voro<BaseV>::get_x(int i, double & x, double & y) const {
  auto p = vh_arr_[i]->point();
  x = p.x();
  y = p.y();
}

template <class BaseV>
template <typename T>
void VM_voro<BaseV>::get_order_para(T &phi, T &theta) const {
  double vx = 0;
  double vy = 0;
  for (int i = 0; i < N_; i++) {
    vx += v_arr_[i].vx;
    vy += v_arr_[i].vy;
  }
  vx /= N_;
  vy /= N_;
  phi = sqrt(vx * vx + vy * vy);
  theta = atan2(vy, vx); 
}