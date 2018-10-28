#include "particle.h"

void Par::update(double noise, double v0, const Vec_2<double>& l) {
  double tmp = sqrt(vx_next * vx_next + vy_next * vy_next);
  double c1 = vx_next / tmp;
  double s1 = vy_next / tmp;
  double c2 = cos(noise);
  double s2 = sin(noise);
  vx = vx_next = c1 * c2 - s1 * s2;
  vy = vy_next = c1 * s2 + c2 * s1;

  x += v0 * vx;
  y += v0 * vy;

  tangle_1(x, l.x);
  tangle_1(y, l.y);

  cell_index = int(x) + int(y) * int(l.x);
}

void Par::update(double noise, double v0, const Vec_2<double>& l,
  const std::vector<Vec_2<double>> &rand_field) {
#ifdef RAND_FIELD
  vx_next += n_neighb * rand_field[cell_index].x;
  vy_next += n_neighb * rand_field[cell_index].y;
  n_neighb = 1;
#endif
  update(noise, v0, l);
}

void Par2::update(double noise, double v0, const Vec_2<double>& l) {
  double tmp = sqrt(vx_next * vx_next + vy_next * vy_next);
  double c1 = vx_next / tmp;
  double s1 = vy_next / tmp;
  double c2 = cos(noise);
  double s2 = sin(noise);
  vx = vx_next = c1 * c2 - s1 * s2;
  vy = vy_next = c1 * s2 + c2 * s1;

  x += v0 * vx;
  y += v0 * vy;

  if (x > l.x) {
    offset.x += 1;
    x -= l.x;
  } else if (x < 0) {
    offset.x -= 1;
    x += l.x;
  }
  if (y > l.y) {
    offset.y += 1;
    y -= l.y;
  } else if (y < 0) {
    offset.y -= 1;
    y += l.y;
  }
  cell_index = int(x) + int(y) * int(l.x);
}

void Par2::update(double noise, double v0, const Vec_2<double>& l, const std::vector<Vec_2<double>>& rand_field) {
#ifdef RAND_FIELD
  vx_next += n_neighb * rand_field[cell_index].x;
  vy_next += n_neighb * rand_field[cell_index].y;
  n_neighb = 1;
#endif
  update(noise, v0, l);
}
