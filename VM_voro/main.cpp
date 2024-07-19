#include "vm_voro.h"

int main(int argc, char* argv[]) {
  //set parameters
 
  double Lx = 100;
  double Ly = 100;
  double eta = 0.2;
  double rho0 = 1;
  double v0 = 0.5;
  int N = int(Lx * Ly * rho0);
  int n_step = 10000;
  int log_dt = 1000;
  int snap_dt = 1000;
  int seed = 3000;
  Ran myran(seed);

  VM_voro<V_scalar> birds(Lx, Ly, N, eta, v0);

  birds.ini_rand(myran, 0);

  for (int i = 1; i < n_step; i++) {
    birds.align();
    birds.stream(myran);
    if (i % 100 == 0) {
      double phi, theta;
      birds.get_order_para(phi, theta);
      std::cout << "t=" << i << "\tphi=" << phi << "\ttheta=" << theta << std::endl;
    }
  }

}
