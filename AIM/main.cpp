#include <iostream>
#include <fstream>
#include <iomanip>
#include "AIM.h"

using namespace std;

int main(int argc, char* argv[]) {
  int Lx = 256;
  int Ly = 100;
  double beta = 1.9;
  double eps = 0.;
  double rho0 = 6;
  double h0 = 0.1;
  int t_half = 100;
  int n_period = 500;
  unsigned long long seed = 5;
  run_osc(Lx, beta, eps, rho0, h0, t_half, n_period, seed, false);
  //run_reverse(Lx, beta, eps, rho0, h0, 2, 200);
  //run_osc(Lx, Ly, beta, eps, rho0, h0, t_half, n_period, seed, true, false);
}