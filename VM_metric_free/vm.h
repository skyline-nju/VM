#ifndef VM_H
#define VM_H
#include <iostream>
#include <vector>
#include <cmath>
#include "particle.h"
#include "rand.h"
#include "comn.h"
#include "cell_list.h"
#include "cmdline.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_2.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_traits_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Vector_2.h>

// Base class for Vicsek model
class VM {
public:
  VM(const cmdline::parser &cmd, Ran &myran);
  virtual ~VM();
 
  // initialize
  virtual void input_data(
    const double *x, const double *y, const double *vx, const double *vy) = 0;
  virtual void ini(const cmdline::parser &cmd, Ran &myran);
  virtual void create_random(Ran &myran);
  virtual void create_defect_pair(
    Ran &myran, double defect_sep, double rot_angle);

  // run
  virtual void align() = 0;
  virtual void align_asym(double alpha) = 0;
  virtual void stream(double dt, Ran &myran) = 0;

  // get data
  virtual void get_x(int i, double &x, double &y) const = 0;
  virtual void get_v(int i, double &vx, double &vy) const = 0;
  virtual void get_v_mean(double &vx_m, double &vy_m) const;
  virtual void get_theta(int i, double &theta) const= 0;
  virtual void get_type(int i, uint32_t& my_type) const = 0;
  int get_num_birds() const { return N; }
  virtual void output_data(double *x, double *y, double *vx, double *vy) const;

protected:
  double v0;
  double Lx;
  double Ly;
  double half_Lx_;
  double half_Ly_;
  int N;
  double eta;
  double extra_torque;
  double *random_torque;
};

// Vicsek mocel for metric case, using cell list to optimize
template <class BaseV>
class VM_metric: public VM {
public:
  typedef ParNode<BaseV> Node;

  VM_metric(const cmdline::parser &cmd, Ran &myran):
    VM(cmd, myran) { ini(cmd, myran); }
  ~VM_metric() { delete[] cell; }
  void input_data(const double *x, const double *y,
                  const double *vx, const double *vy);
  void align();
  void align_asym(double alpha){};
  void stream(double dt, Ran &myran);
  void get_x(int i, double &x, double &y) const { x = par_arr[i].x; y = par_arr[i].y; }
  void get_v(int i, double &vx, double &vy) const { par_arr[i].get_v(vx, vy); }
  void get_theta(int i, double &Theta) const { par_arr[i].get_theta(Theta); }
  void get_type(int i, uint32_t& my_type) const { my_type = 0; }

  std::vector<Node> par_arr;
  Cell<Node> *cell;
};

template<class BaseV>
inline void VM_metric<BaseV>::input_data(const double * x, const double * y,
                                         const double * vx, const double * vy) {
  par_arr.reserve(N);
  for (int i = 0; i < N; i++) {
    par_arr.emplace_back(x[i], y[i], vx[i], vy[i]);
  }
  Cell<Node>::mx = int(Lx);
  Cell<Node>::my = int(Ly);
  Cell<Node>::mm = Cell<Node>::mx * Cell<Node>::my;
  Cell<Node>::l0 = 1;
  cell = new Cell<Node>[Cell<Node>::mm];
}

template <class BaseV>
inline void VM_metric<BaseV>::align() {
  Cell<Node>::refresh(cell, par_arr);
  Cell<Node>::all_pairs(cell, Lx, Ly);
}

template <class BaseV>
void VM_metric<BaseV>::stream(double dt, Ran &myran) {
  for (auto it = par_arr.begin(); it != par_arr.end(); ++it) {
    double torque = extra_torque;
    if (random_torque) {
      torque += random_torque[(*it).cell_idx];
    }
    (*it).move(eta, torque, myran, v0, Lx, Ly);
  }
}

// Vicsek model for metric-free case, using CGAL.
template <class BaseV>
class VM_metric_free : public VM {
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
  typedef std::pair<Point, unsigned int> Pair_P_I;

  VM_metric_free(const cmdline::parser &cmd, Ran &myran);
  void input_data(const double *x, const double *y,
                  const double *vx, const double *vy);
  void align();
  void align_asym(double alpha);
  void stream(double dt, Ran &myran);
  void get_x(int i, double &x, double &y) const;
  void get_v(int i, double &vx, double &vy) const;
  void get_theta(int i, double &theta)const;
  void get_dR(int i, int j, double &dx, double &dy) const;
  void get_type(int i, uint32_t& type_i) const { v_arr[i].get_type(type_i); }

  std::vector<BaseV> v_arr;
  std::vector<Pair_P_I> x_arr;
  PDT *DT;
};

template<class BaseV>
VM_metric_free<BaseV>::VM_metric_free(const cmdline::parser & cmd, Ran & myran):
                                       VM(cmd, myran) {
  // initialize triangulation
  PDT::Iso_rectangle domain(0, 0, Lx, Ly);
  DT = new PDT(domain);
  
  // initialize v_arr and x_arr
  ini(cmd, myran);
  if (cmd.exist("dis_frac")) {
    int n_dis = int(N * cmd.get<double>("dis_frac"));
    for (int i = 0;  i < n_dis; i++) {
      v_arr[i].set_type(1);
    }
    std::cout << "there are " << n_dis << " dissenters out of " << N << " birds" << std::endl;
  }
}

template<class BaseV>
void VM_metric_free<BaseV>::input_data(const double * x, const double * y,
                                       const double * vx, const double * vy) {
  v_arr.reserve(N);
  x_arr.reserve(N);
  for (int i = 0; i < N; i++) {
    v_arr.emplace_back(vx[i], vy[i]);
    x_arr.push_back(std::make_pair(Point(x[i], y[i]), i));
  }
}

template<class BaseV>
void VM_metric_free<BaseV>::align() {
  // Bug: fail to insert particles when Lx > Ly
  DT->insert(x_arr.begin(), x_arr.end());
  // Current solution is avoiding Lx > Ly
  
  for (auto fit = DT->periodic_triangles_begin(PDT::UNIQUE);
    fit != DT->periodic_triangles_end(PDT::UNIQUE); ++fit) {
    unsigned int idx0 = fit.get_face()->vertex(0)->info();
    unsigned int idx1 = fit.get_face()->vertex(1)->info();
    unsigned int idx2 = fit.get_face()->vertex(2)->info();
    if (idx0 < idx1) {
      v_arr[idx0].collide(&v_arr[idx1]);
    }
    if (idx1 < idx2) {
      v_arr[idx1].collide(&v_arr[idx2]);
    }
    if (idx2 < idx0) {
      v_arr[idx2].collide(&v_arr[idx0]);
    }
  }
  DT->clear();
}

template<class BaseV>
void VM_metric_free<BaseV>::align_asym(double alpha) {
  // Bug: fail to insert particles when Lx > Ly
  DT->insert(x_arr.begin(), x_arr.end());
  // Current solution is avoiding Lx > Ly
  // Besides, the inserting speed is very slow for unequal Lx and Ly
  // such that for now it seems we need run with Lx=Ly
  
  for (auto fit = DT->periodic_triangles_begin(PDT::UNIQUE);
    fit != DT->periodic_triangles_end(PDT::UNIQUE); ++fit) {
    unsigned int idx0 = fit.get_face()->vertex(0)->info();
    unsigned int idx1 = fit.get_face()->vertex(1)->info();
    unsigned int idx2 = fit.get_face()->vertex(2)->info();
    double dx, dy;
    if (idx0 < idx1) {
      get_dR(idx0, idx1, dx, dy);
      v_arr[idx0].collide_asym(&v_arr[idx1], dx, dy, alpha);
    }
    if (idx1 < idx2) {
      get_dR(idx1, idx2, dx, dy);
      v_arr[idx1].collide_asym(&v_arr[idx2], dx, dy, alpha);
    }
    if (idx2 < idx0) {
      get_dR(idx2, idx0, dx, dy);
      v_arr[idx2].collide_asym(&v_arr[idx0], dx, dy, alpha);
    }
  }
  DT->clear();
}
template<class BaseV>
void VM_metric_free<BaseV>::get_dR(int i, int j, double &dx, double &dy) const {
  dx = (x_arr[j].first)[0] - (x_arr[i].first)[0];
  dy = (x_arr[j].first)[1] - (x_arr[i].first)[1];
  if (dx > half_Lx_) {
    dx -= Lx;
  } else if (dx < -half_Lx_) {
    dx += Lx;
  }
  if (dy > half_Ly_) {
    dy -= Ly;
  } else if (dy < -half_Ly_) {
    dy += Ly;
  }
}

template<class BaseV>
void VM_metric_free<BaseV>::stream(double dt, Ran & myran) {
  for (int i = 0; i < N; i++) {
    double torque = extra_torque;
    if (random_torque) {
      double x, y;
      get_x(i, x, y);
      int col = int(x);
      int row = int(y);
      int idx = col + row * int(Lx);
      torque += random_torque[idx];
    }
    v_arr[i].update_v(eta, torque, myran);
    double x = (x_arr[i].first)[0];
    double y = (x_arr[i].first)[1];
    double x0 = x;
    double y0 = y;
    v_arr[i].update_x(x, y, v0, Lx, Ly);
    x_arr[i].first += Vec_2(x - x0, y - y0);
  }
}

template<class BaseV>
inline void VM_metric_free<BaseV>::get_x(int i, double & x, double & y) const{
  x = x_arr[i].first.x();
  y = x_arr[i].first.y();
}

template<class BaseV>
inline void VM_metric_free<BaseV>::get_v(int i, double & vx, double & vy) const{
  v_arr[i].get_v(vx, vy);
}

template<class BaseV>
inline void VM_metric_free<BaseV>::get_theta(int i, double & Theta) const{
  v_arr[i].get_theta(Theta);
}

#endif