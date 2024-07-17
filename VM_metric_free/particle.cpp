#include "particle.h"

void V_scalar::set_v(double theta) {
  vx = cos(theta);
  vy = sin(theta);
  vx_next = vx;
  vy_next = vy;
}

void V_scalar::update_v(double eta, double torque, Ran &myran) {
  double tmp = std::sqrt(vx_next * vx_next + vy_next * vy_next);
  double c1 = vx_next / tmp;
  double s1 = vy_next / tmp;
  double noise;
  if (eta > 0) {
    noise = (myran.doub() - 0.5) * eta * 2 * PI + torque;
  } else {
    noise = torque;
  }
  double c2 = std::cos(noise);
  double s2 = std::sin(noise);
  vx = vx_next = c1 * c2 - s1 * s2;
  vy = vy_next = c1 * s2 + c2 * s1;
}

void V_scalar::update_x(double & _x, double & _y, double v0, double Lx, double Ly) {
  _x += v0 * vx;
  _y += v0 * vy;
  if (_x < 0) {
    _x += Lx;
  } else if (_x >= Lx) {
    _x -= Lx;
  }
  if (_y < 0) {
    _y += Ly;
  } else if (_y >= Ly) {
    _y -= Ly;
  }
}


void V_scalar_aligner_dissenter::set_v(double theta) {
  if (is_aligner) {
      vx = cos(theta);
      vy = sin(theta);
      vx_next = vx;
      vy_next = vy;
  }
}

void V_vectorial::update_v(double eta, double torque, Ran &myran) {
  double noise_x, noise_y;
  myran.circle_point_picking(noise_x, noise_y);
  double tmp = eta * n_neighbor;
  vx_next += noise_x * tmp;
  vy_next += noise_y * tmp;
  tmp = 1 / std::sqrt(vx_next * vx_next + vy_next * vy_next);
  vx_next *= tmp;
  vy_next *= tmp;
  if (torque != 0) {
    double c1 = std::cos(torque);
    double s1 = std::sin(torque);
    vx = c1 * vx_next - s1 * vy_next;
    vy = c1 * vy_next - s1 * vx_next;
    vy_next = vx;
    vy_next = vy;
  } else {
    vx = vx_next;
    vy = vy_next;
  }
  n_neighbor = 1;
}

double V_conti::h;
double V_conti::sqrt_24_Dr_h;

void V_conti::update_x(double & _x, double & _y, double v0, double Lx, double Ly) {
  double vx = v0 * std::cos(theta);
  double vy = v0 * std::sin(theta);
  _x += vx * h;
  _y += vy * h;
  if (_x < 0) {
    _x += Lx;
  } else if (_x >= Lx) {
    _x -= Lx;
  }
  if (_y < 0) {
    _y += Ly;
  } else if (_y >= Ly) {
    _y -= Ly;
  }
}
