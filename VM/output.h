#ifndef OUTPUT_H
#define OUTPUT_H
#include <ctime>
#include "grid.h"
#ifdef _MSC_VER
#include <io.h>
#else
#include <unistd.h>
#endif


void output_ini(double eta, double epsilon, unsigned long long seed, int nStep);

void output(const Node *bird, int step);

void mkdir(const char *folder);

#endif
