#include <iostream>
#include <fstream>
#include <iomanip>
#include "ACMq4.h"

using namespace std;

int main(int argc, char* argv[]) {
  int Lx = 64;
  int Ly = 256;
  double beta = 4.;
  double eps = 0.9;
  double rho0 = 6;
  double D = 0.1;
  int n_step = 10000000;
  int dn_out = 5000;
  unsigned long long seed = 3001;


  run(Lx, Ly, rho0, beta, eps, D, n_step, dn_out, seed);
}