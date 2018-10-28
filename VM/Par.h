#ifndef NODE_H
#define NODE_H
#include "rand.h"
#include "comn.h"


struct Par
{
  Par();
  double rr(Par *node);
  double rr(Par *node, double a, double b);
  void addV(Par *node);
  void align(Par *node);
  void align(Par *node, double a, double b);
  void move(double noise);
  void move(double noise, Ran *myran);

  static Par *ini_rand(Ran *myran);
  static Par *ini_move_left(Ran *myran);
  static Par *ini_from_snap(double Lx0, double Ly0,
                             const std::vector<float> &x0,
                             const std::vector<float> &y0,
                             const std::vector<float> &theta0);
  double x, y, vx, vx_next, vy, vy_next;
  int cell_idx;
  Par* next;

  static double Lx;
  static double Ly;
  static double v0;
  static double rho_0;
  static int N;
};

inline Par::Par() {
  x = y = vx = vx_next = vy = vy_next;
  cell_idx = 0;
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

#endif