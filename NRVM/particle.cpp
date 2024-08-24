#include "particle.h"
#include <cmath>

Par::Par(double x0, double y0, double vx0, double vy0) {
  x = x0;
  y = y0;
  vx = vx_next = vx0;
  vy = vy_next = vy0;
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


