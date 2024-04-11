#include <iostream>
#include <fstream>
#include <iomanip>
#include "ACMq4_QS.h"

using namespace std;

int main(int argc, char* argv[]) {
  int Lx = 128;
  int Ly = 128;
  double beta = 3.2;
  double eps = 1;
  double rho0 = 10;
  double eta = -0.8;
  double D = 0.1;
  int n_step = 10000000;
  int dn_out = 5000;
  unsigned long long seed = 3001;

  double rho_thresh = 10;


  //run(Lx, Ly, rho0, beta, eps, D, n_step, dn_out, seed);
  run_QS(Lx, Ly, rho0, beta, eps, D, rho_thresh, eta, n_step, dn_out, seed);

}