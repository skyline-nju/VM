#ifndef OUTPUT_H
#define OUTPUT_H
#include <ctime>
#include "node.h"

void output_ini(double eta, double epsilon, unsigned long long seed, int nStep, int nCell);

void output(const Node *bird, int step);


#endif
