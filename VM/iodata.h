#ifndef IODATA_H
#define IODATA_H
#include "node.h"

void ini_output(double eta, double epsilon, unsigned long long seed, int nStep, int nCell, bool flag_snap_one, int snap_interval);

void output(const Node *bird, int step);
#endif
