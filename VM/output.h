#ifndef OUTPUT_H
#define OUTPUT_H
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <ctime>
#ifdef _MSC_VER
#include <io.h>
#else
#include <unistd.h>
#endif
#include "grid.h"

void output_ini(double eta, double epsilon, int seed, int Nstep);

void output(const Node *bird, int step);

void mkdir(char *folder);
#endif
