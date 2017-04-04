#ifndef RUN_H
#define RUN_H
#include "grid.h"
#include "cmdline.h"
#include "iodata.h"

void ini_birds(Node **p_bird, Ran *myran, const cmdline::parser &cmd);

void ini_rand_torques(double **disorder, int n, double epsilon, unsigned long long seed);

void run(Node *bird, Grid *cell, Ran*myran, int nStep, double eta);

void run(Node *bird, Grid *cell, Ran*myran, int nStep, double eta, double epsilon, const double *disorder);

#endif
