#include <iostream>
#include "comn.h"
#include "rand.h"
#include "node.h"
#include "particle.h"
#include "cellList.h"
#include "configure.h"
#include "run.h"

int main(int argc, char* argv[]) {
#ifdef _MSC_VER
  double eta = 0.18;
  double eps = 0;
  int L = 32;
  int seed = 1237;
  int n_step = 10000000;
#else
  double eta = atof(argv[1]);
  double eps = atof(argv[2]);
  double L = atof(argv[3]);
  int n_step = atoi(argv[4]);
  int seed = atoi(argv[5]);
#endif
  double rho_0 = 1.;
  double v0 = 0.5;
  double eta2PI = eta * 2 * PI;

  int n_par = int(L * L * rho_0);
  Ran myran(seed);
  Vec_2<double> gl_l(L, L);
  Vec_2<double> origin(0, 0);
  Vec_2<double> half_l(L * 0.5, L * 0.5);

#ifdef CAL_MSD
  typedef UniNode<Par2> par_t;
#else
  typedef UniNode<Par> par_t;
#endif

  std::vector<par_t> p_arr;
  p_arr.reserve(n_par);
  for (int i = 0; i < n_par; i++) {
    p_arr.emplace_back(myran, gl_l, origin);
  }

  CellListNode_2<par_t> cl(gl_l, 1., origin);
  cl.create(p_arr);

  auto pair_force = [&gl_l, &half_l](par_t *pi, par_t *pj) {
    if (pi->rr(pj, gl_l, half_l) < 1) {
      pi->addV(pj);
    }
  };

#ifdef RAND_FIELD
  std::vector<Vec_2<double>> rand_field;
  ini_rand_field(rand_field, eps, n_par, myran);
  auto move = [n_par, &myran, eta2PI, v0, &gl_l, &p_arr, &rand_field, &cl]() {
    for (int i = 0; i < n_par; i++) {
      const double noise = (myran.doub() - 0.5) * eta2PI;
      p_arr[i].update(noise, v0, gl_l, rand_field);
    }
    cl.recreate(p_arr);
  };
#else
  std::vector<double> rand_torque;
  ini_rand_torque(rand_torque, eps, n_par, myran);
  auto move = [n_par, &myran, eta2PI, v0, &gl_l, &p_arr, &rand_torque, &cl]() {
    for (int i = 0; i < n_par; i++) {
      const double noise = (myran.doub() - 0.5) * eta2PI + rand_torque[p_arr[i].cell_index];
      p_arr[i].update(noise, v0, gl_l);
    }
    cl.recreate(p_arr);
  };
#endif

#ifdef CAL_MSD
  run_MSD(eta, eps, gl_l, seed, n_step, pair_force, move, p_arr, cl);
#else
  run(eta, eps, gl_l, seed, n_step, pair_force, move, p_arr, cl);
#endif
}