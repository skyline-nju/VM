#pragma once

#include "config.h"
#include "vm.h"

#ifdef MOVE_POINTS_ON

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
  typedef Tds::Vertex_handle Vertex_handle;

  VM_metric_free(const cmdline::parser &cmd, Ran &myran);
  ~VM_metric_free() { delete DT; }
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
  void set_v(double theta);

  std::vector<BaseV> v_arr;
  std::vector<Vertex_handle> vh_arr;
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
  vh_arr.reserve(N);
  for (int i = 0; i < N; i++) {
    v_arr.emplace_back(vx[i], vy[i]);
    auto vh = DT->insert(Point(x[i], y[i]));
    vh->info() = i;
    vh_arr.push_back(vh);
  }
}


template <class BaseV>
void VM_metric_free<BaseV>::set_v(double theta) {
  for (auto& v: v_arr) {
    v.set_v(theta);
  }
}
template<class BaseV>
void VM_metric_free<BaseV>::align() {  
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
}

template<class BaseV>
void VM_metric_free<BaseV>::align_asym(double alpha) {
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
}

template<class BaseV>
void VM_metric_free<BaseV>::get_dR(int i, int j, double &dx, double &dy) const {
  double xi, yi, xj, yj;
  get_x(i, xi, yi);
  get_x(j, xj, yj);
  dx = xj - xi;
  dy = yj - yi;
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
    v_arr[i].update_v(eta, extra_torque, myran);
    double x, y;
    get_x(i, x, y);
    v_arr[i].update_x(x, y, v0, Lx, Ly);
    Point p_new(x, y);
    auto vh_old = vh_arr[i];
    vh_arr[i] = DT->insert(p_new, vh_old->face());
    vh_arr[i]->info() = vh_old->info();
    DT->remove(vh_old);
  }

}

template<class BaseV>
inline void VM_metric_free<BaseV>::get_x(int i, double & x, double & y) const{
  auto p = vh_arr[i]->point();
  x = p.x();
  y = p.y();
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