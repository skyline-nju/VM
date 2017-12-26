#ifndef METRIC_FREE_VM_H
#define METRIC_FREE_VM_H
#include <fstream>
#include <ctime>
#include <chrono>
#include "bird.h"

void ini_rand(double Lx, double Ly, double eta,
              unsigned long long seed, int nStep);

void run();

void run(int n);

void run(int n, double &phi);

void run_test(int n);

#endif