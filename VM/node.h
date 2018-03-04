#ifndef NODE_H
#define NODE_H
#include "rand.h"
#include "comn.h"


struct Node
{
  Node();
  double rr(Node *node);
  double rr(Node *node, double a, double b);
  void addV(Node *node);
  void align(Node *node);
  void align(Node *node, double a, double b);
  void move(double noise);
  void move(double noise, Ran *myran);

  static Node *ini_rand(Ran *myran);
  static Node *ini_move_left(Ran *myran);
  static Node *ini_from_snap(double Lx0, double Ly0,
                             const std::vector<float> &x0,
                             const std::vector<float> &y0,
                             const std::vector<float> &theta0);
  double x, y, vx, vx0, vy, vy0;
  int cell_idx;
  Node* next;

  static double Lx;
  static double Ly;
  static double v0;
  static double rho_0;
  static int N;
};

inline Node::Node() {
  x = y = vx = vx0 = vy = vy0;
  cell_idx = 0;
  next = nullptr;
}

inline double Node::rr(Node *node) {
  double dx = node->x - x;
  double dy = node->y - y;
  return dx*dx + dy*dy;
}

inline double Node::rr(Node *node, double a, double b) {
  double dx = node->x - x + a;
  double dy = node->y - y + b;
  return dx*dx + dy*dy;
}

inline void Node::addV(Node *node) {
  vx += node->vx0;
  vy += node->vy0;
  node->vx += vx0;
  node->vy += vy0;
}

#endif