#ifndef PARTICLE_H
#define PARTICLE_H
#include <cmath>
#include "rand.h"
#include "comn.h"

// velocity of particles that move in discrete time steps and suffer the 
// scalar noise.
struct V_scalar {
  V_scalar() {}
  V_scalar(double _vx, double _vy) :
    vx(_vx), vy(_vy), vx_next(_vx), vy_next(_vy) {
  }

  template <class T>
  void collide(T *other);
  void update_v(double eta, double torque, Ran &myran);
  void update_x(double &_x, double &_y, double v0, double Lx, double Ly);
  void get_v(double &VX, double &VY) const { VX = vx; VY = vy; }
  void get_theta(double &Theta) const { Theta = std::atan2(vy, vx); }

  double vx;
  double vy;
  double vx_next;
  double vy_next;
};

template <class T>
inline void V_scalar::collide(T *other) {
  vx_next += other->vx;
  vy_next += other->vy;
  other->vx_next += vx;
  other->vy_next += vy;
}

// velocity of particles that move in discrete time steps and suffer the 
// vectorial noise.
struct V_vectorial : public V_scalar {
  V_vectorial() {}
  V_vectorial(double _vx, double _vy) :
    V_scalar(_vx, _vy), n_neighbor(1) {}

  template <class T>
  void collide(T *other);
  void update_v(double eta, double torque, Ran &myran);

  int n_neighbor;
};

template <class T>
inline void V_vectorial::collide(T *other) {
  vx_next += other->vx;
  vy_next += other->vy;
  other->vx_next += vx;
  other->vy_next += vy;
  n_neighbor++;
  other->n_neighbor++;
}

// orientation of particles that move in continuos time steps.
struct V_conti {
  V_conti() {}
  V_conti(double vx, double vy):
    theta(std::atan2(vy, vx)), theta_dot(0) {}

  template <class T>
  void collide(T *other);

  void update_v(double eta, double torque, Ran &myran) {
    theta += (theta_dot + torque) * h + (myran.doub() - 0.5) * sqrt_24_Dr_h;
    theta_dot = 0;
  }
  void update_x(double &_x, double &_y, double v0, double Lx, double Ly);

  void get_v(double &VX, double &VY) const {
    VX = std::cos(theta); VY = std::sin(theta);
  }

  void get_theta(double &Theta) const { Theta = theta; }

  double theta;
  double theta_dot;
  static double h;               // time interval to integrate
  static double sqrt_24_Dr_h;    // sqrt(2 * Dr * 12 * h)
};

template <class T>
inline void V_conti::collide(T *other) {
  double omega = std::sin(other->theta - theta);
  theta_dot += omega;
  other->theta_dot -= omega;
}

// Linknode for particle moving with metric-interaction
template <class BaseV>
class ParNode :public BaseV {
public:
  ParNode() : BaseV() {}
  ParNode(double _x, double _y, double _vx, double _vy) :
    BaseV(_vx, _vy), x(_x), y(_y), next(nullptr) {
  }

  bool within_range(const ParNode<BaseV> *other) const;
  bool within_range(const ParNode<BaseV> *other, double a, double b) const;
  void interact(ParNode *other);
  void interact(ParNode *other, double a, double b);
  void move(double eta, double torque, Ran &myran, double v0, double Lx, double Ly);

  double x;
  double y;
  int cell_idx;
  ParNode<BaseV> *next;
};

template <class BaseV>
bool ParNode<BaseV>::within_range(const ParNode<BaseV> *other) const {
  double dx = other->x - x;
  double dy = other->y - y;
  return (dx * dx + dy * dy < 1);
}

template <class BaseV>
bool ParNode<BaseV>::within_range(const ParNode<BaseV> *other, double a, double b) const {
  double dx = other->x - x + a;
  double dy = other->y - y + b;
  return (dx * dx + dy * dy < 1);
}

template <class BaseV>
inline void ParNode<BaseV>::interact(ParNode *other) {
  if (within_range(other))
    this->collide(other);
}

template <class BaseV>
inline void ParNode<BaseV>::interact(ParNode *other, double a, double b) {
  if (within_range(other, a, b))
    this->collide(other);
}

template <class BaseV>
inline void ParNode<BaseV>::move(double eta, double torque, Ran &myran,
                                 double v0, double Lx, double Ly) {
  this->update_v(eta, torque, myran);
  this->update_x(x, y, v0, Lx, Ly);
}

#endif