#include <iostream>
#include <fstream>
#include <iomanip>
#include "nematic_lattice_gas.h"

using namespace std;

int main(int argc, char* argv[]) {
  int Lx = 256;
  int Ly = 128;
  double beta = 2;
  double eps = 0.9;
  double rho0 = 1.85;
  double h0 = 0.1;
  int n_period = 100000;
  int output_interval = 10000;
  unsigned long long seed = 7;
  double D = 1;

  run(Lx, Ly, beta, eps, rho0, D, seed, n_period, output_interval);
}