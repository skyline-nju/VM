#ifndef NODE_H
#define NODE_H
#include "rand.h"
#include "comn.h"
#include <vector>
#include <iostream>

template <typename TPar>
class Node: public TPar {
public:
  Node(): TPar(), next(nullptr) {}
  Node(double x0, double y0, double vx0, double vy0):
    TPar(x0, y0, vx0, vy0), next(nullptr) {}

  Node<TPar> *next;
};

class Par {
public:
  Par();

  Par(double x0, double y0, double vx0, double vy0);


  template <typename TPar>
  double rr(TPar *node);

  template <typename TPar>
  double rr(TPar *node, double a, double b);


  template <typename TPar>
  void addV(TPar *node);

  template <typename TPar>
  void align(TPar *node);

  template <typename TPar>
  void align(TPar *node, double a, double b);

  void move(double noise, double Lx, double Ly, double v0);

  double get_theta() const {
    return atan2(vy, vx);
  }

  template<typename TRan>
  void move(double eta, TRan& myran, double Lx, double Ly, double v0);


  void set_type (int type_id0) {}

  int get_type() {return 0;}

  double x, y, vx, vx_next, vy, vy_next;
};

inline Par::Par() {
  x = y = vx = vx_next = vy = vy_next = 0;
}

template<typename TPar>
inline double Par::rr(TPar *node) {
  double dx = node->x - x;
  double dy = node->y - y;
  return dx*dx + dy*dy;
}

template<typename TPar>
inline double Par::rr(TPar *node, double a, double b) {
  double dx = node->x - x + a;
  double dy = node->y - y + b;
  return dx*dx + dy*dy;
}

template<typename TPar>
inline void Par::addV(TPar *node) {
  vx += node->vx_next;
  vy += node->vy_next;
  node->vx += vx_next;
  node->vy += vy_next;
}

template<typename TPar>
void Par::align(TPar *node) {
  if (rr(node) < 1)
    addV(node);
}

template<typename TPar>
void Par::align(TPar *node, double a, double b) {
  if (rr(node, a, b) < 1)
    addV(node);
}

template<typename TRan>
void Par::move(double eta, TRan& myran, double Lx, double Ly, double v0) {
  double noise = (myran.doub() - 0.5) * 2 * PI * eta;
  move(noise, Lx, Ly, v0);
}

template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}


class Par_w_type: public Par {
public:
  Par_w_type(double x0, double y0, double vx0, double vy0): Par(x0, y0, vx0, vy0), type_id(0) {}

  void set_type (int type_id0) {type_id = type_id0;}

  int get_type() {return type_id;}

  template <typename TPar>
  void addV_NR(TPar *node);

  template<typename TPar>
  void align(TPar *node);

  template<typename TPar>
  void align(TPar *node, double a, double b);

  int type_id;
};

template <typename TPar>
void Par_w_type::addV_NR(TPar *node) {
  if (type_id == 0) {
    vx += node->vx_next;
    vy += node->vy_next;
  }
  if (node->type_id == 0) {
    node->vx += vx_next;
    node->vy += vy_next;
  }
}

template<typename TPar>
void Par_w_type::align(TPar *node) {
  if (rr(node) < 1)
    addV_NR(node);
}

template<typename TPar>
void Par_w_type::align(TPar *node, double a, double b) {
  if (rr(node, a, b) < 1)
    addV_NR(node);
}


#endif


