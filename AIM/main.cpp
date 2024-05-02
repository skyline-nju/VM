#include <iostream>
#include <fstream>
#include <iomanip>
#include "AIM.h"

using namespace std;

int main(int argc, char* argv[]) {
  int Lx = 512;
  int Ly = 32;
  double beta = 1.8;
  double eps = 0.9;
  double rho0 = 3;
  double D = 1;
  double alpha = 0;
  int n_step = 1000000;
  int dn_out = 1000;
  unsigned long long seed = 4002;
  std::string ini_condi = "ordered";

  run(Lx, Ly, rho0, beta, eps, D, n_step, dn_out, seed, alpha, ini_condi);
}