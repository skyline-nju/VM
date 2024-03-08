#include <iostream>
#include <fstream>
#include <iomanip>
#include "AIM.h"

using namespace std;

int main(int argc, char* argv[]) {
  int Lx = 1024;
  int Ly = 64;
  double beta = 5;
  double eps = 1;
  double rho0 = 5;
  double D = 1;
  double alpha = 1;
  int n_step = 10000000;
  int dn_out = 5000;
  unsigned long long seed = 1001;


  run(Lx, Ly, rho0, beta, eps, D, n_step, dn_out, seed, alpha);
}