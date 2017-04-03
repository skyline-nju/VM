#ifndef RUN_H
#define RUN_H
#include "grid.h"

void ini_rand_torques(double **disorder, int n, double epsilon, Ran *myran);

void run(Node *bird, Grid *cell, Ran*myran, int nStep, double eta);

void run(Node *bird, Grid *cell, Ran*myran, int nStep, double eta, double epsilon, const double *disorder);

#endif
