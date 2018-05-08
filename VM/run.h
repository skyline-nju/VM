#ifndef RUN_H
#define RUN_H
#include "grid.h"
#include "cmdline.h"
#include "exporter.h"

void ini_output(const cmdline::parser &cmd);

void ini_birds(Node **p_bird, int &n_par, Ran *myran,
               const cmdline::parser &cmd);

void ini_rand_torques(
  double **disorder, int n, double epsilon, unsigned long long seed);

void run(Node *bird, int n_bird, Grid *cell, Ran *myran, int nStep, double eta,
         const double *disorder, bool vicsekShake=false);

void finish_simulation();
#endif
