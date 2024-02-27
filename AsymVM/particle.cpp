#include "particle.h"
#include <cmath>

double Par::alpha = 1.;

Par::Par(double x0, double y0, double theta0) {
  x = x0;
  y = y0;
  vx = vx_next = cos(theta0);
  vy = vy_next = sin(theta0);
  next = nullptr;
}

void Par::align(Par *node) {
  if (rr(node) < 1)
    addV(node);
}

void Par::align(Par *node, double a, double b) {
  if (rr(node, a, b) < 1)
    addV(node);
}

void Par::asym_align(Par* node) {
  double dx = node->x - x;
  double dy = node->y - y;
  double rr = dx * dx + dy * dy;

  if (rr < 1) {
    double w1 = 0.5 * (1 + alpha * sgn(dx * vx_next + dy * vy_next));
    double w2 = 0.5 * (1 + alpha * sgn(-dx * node->vx_next - dy * node->vy_next));
    addV(node, w1, w2);
  }
}

void Par::asym_align(Par* node, double a, double b) {
  double dx = node->x - x + a;
  double dy = node->y - y + b;
  double rr = dx * dx + dy * dy;

  if (rr < 1) {
    double w1 = 0.5 * (1 + alpha * sgn(dx * vx_next + dy * vy_next));
    double w2 = 0.5 * (1 + alpha * sgn(-dx * node->vx_next - dy * node->vy_next));
    addV(node, w1, w2);
  }
}


void Par::move(double noise, double Lx, double Ly, double v0) {
  double tmp = sqrt(vx*vx + vy*vy);
  double c1 = vx / tmp;
  double s1 = vy / tmp;
  double c2 = cos(noise);
  double s2 = sin(noise);
  vx = vx_next = c1 * c2 - s1 * s2;
  vy = vy_next = c1 * s2 + c2 * s1;

  x += v0*vx;
  if (x >= Lx)
    x -= Lx;
  else if (x < 0)
    x += Lx;
  y += v0*vy;
  if (y >= Ly)
    y -= Ly;
  else if (y < 0)
    y += Ly;

}
