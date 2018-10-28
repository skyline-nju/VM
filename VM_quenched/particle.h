#pragma once
#include "configure.h"
#include "vect.h"
#include <cstdlib>
#include "comn.h"
class Par {
public:
  Par() {}

  template<typename TRan>
  Par(TRan &myran, const Vec_2<double> &l, const Vec_2<double> &origin);

  double rr(Par* node);

  double rr(Par* node, double a, double b);

  double rr(Par* node, const Vec_2<double> &l, const Vec_2<double> &half_l);

  void addV(Par* p2);

  void update(double noise, double v0, const Vec_2<double> &l);

  void update(double noise, double v0, const Vec_2<double> &l,
              const std::vector<Vec_2<double>> &rand_field);

  double x = 0;
  double y = 0;
  double vx = 0;
  double vy = 0;
  double vx_next = 0;
  double vy_next = 0;
  int cell_index = 0;
#ifdef RAND_FIELD
  int n_neighb  = 1;
#endif
};

template <typename TRan>
Par::Par(TRan& myran, const Vec_2<double>& l, const Vec_2<double>& origin) {
  x = myran.doub() * l.x + origin.x;
  y = myran.doub() * l.y + origin.y;
  circle_point_picking(vx, vy, myran);
  vx_next = vx;
  vy_next = vy;
  cell_index = int(x) + int(y) * int(l.x);
}

inline double Par::rr(Par* node) {
  double dx = node->x - x;
  double dy = node->y - y;
  return dx * dx + dy * dy;
}

inline double Par::rr(Par* node, double a, double b) {
  double dx = node->x - x + a;
  double dy = node->y - y + b;
  return dx * dx + dy * dy;
}

inline double Par::rr(Par* node, const Vec_2<double> &l, const Vec_2<double> &half_l) {
  double dx = node->x - x;
  untangle_1(dx, l.x, half_l.x);
  double dy = node->y - y;
  untangle_1(dy, l.y, half_l.y);
  return dx * dx + dy * dy;
}

inline void Par::addV(Par* p2) {
  vx_next += p2->vx;
  vy_next += p2->vy;
  p2->vx_next += vx;
  p2->vy_next += vy;
#ifdef RAND_FIELD
  n_neighb++;
  p2->n_neighb++;
#endif
}

class Par2: public Par {
public:
  Par2() {}

  template<typename TRan>
  Par2(TRan &myran, const Vec_2<double> &l, const Vec_2<double> &origin): Par(myran, l, origin) {}

  void update(double noise, double v0, const Vec_2<double> &l);

  void update(double noise, double v0, const Vec_2<double> &l,
              const std::vector<Vec_2<double>> &rand_field);

  Vec_2<int> offset{};
};