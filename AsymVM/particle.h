#ifndef NODE_H
#define NODE_H
#include "rand.h"
#include "comn.h"
#include <vector>


struct Par {
  Par();

  Par(double x0, double y0, double theta0);

  template<typename TRan>
  Par(TRan& myran, double Lx, double Ly);

  template<typename TRan>
  Par(TRan& myran, double Lx, double Ly, double theta0);

  double rr(Par *node);
  double rr(Par *node, double a, double b);

  void addV(Par *node);
  void addV(Par* node, double w1, double w2);

  void align(Par *node);
  void align(Par *node, double a, double b);

  void asym_align(Par* node);
  void asym_align(Par* node, double a, double b);

  void move(double noise, double Lx, double Ly, double v0);

  double get_theta() const {
    return atan2(vy, vx);
  }

  template<typename TRan>
  void move(double eta, TRan& myran, double Lx, double Ly, double v0);

  double x, y, vx, vx_next, vy, vy_next;
  Par* next;

  static double alpha;
};

inline Par::Par() {
  x = y = vx = vx_next = vy = vy_next = 0;
  next = nullptr;
}

template<typename TRan>
Par::Par(TRan& myran, double Lx, double Ly) {
  x = myran.doub() * Lx;
  y = myran.doub() * Ly;
  double theta = myran.doub() * 2 * PI;
  vx = cos(theta);
  vy = sin(theta);
  vx_next = vx;
  vy_next = vy;
  next = nullptr;
}

template<typename TRan>
Par::Par(TRan& myran, double Lx, double Ly, double theta0) {
  x = myran.doub() * Lx;
  y = myran.doub() * Ly;
  vx = cos(theta0);
  vy = sin(theta0);
  vx_next = vx;
  vy_next = vy;
  next = nullptr;
}

inline double Par::rr(Par *node) {
  double dx = node->x - x;
  double dy = node->y - y;
  return dx*dx + dy*dy;
}

inline double Par::rr(Par *node, double a, double b) {
  double dx = node->x - x + a;
  double dy = node->y - y + b;
  return dx*dx + dy*dy;
}

inline void Par::addV(Par *node) {
  vx += node->vx_next;
  vy += node->vy_next;
  node->vx += vx_next;
  node->vy += vy_next;
}

inline void Par::addV(Par* node, double w1, double w2) {
  vx += node->vx_next * w1;
  vy += node->vy_next * w1;
  node->vx += vx_next * w2;
  node->vy += vy_next * w2;
}

template<typename TRan>
void Par::move(double eta, TRan& myran, double Lx, double Ly, double v0) {
  double noise = (myran.doub() - 0.5) * 2 * PI * eta;
  move(noise, Lx, Ly, v0);
}



template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

#endif


