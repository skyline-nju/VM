#include "metric_free_vm.h"

double Lx;
double Ly;
double eta;
int nStep;
double rho0 = 1;
double v0 = 0.5;

Ran *myran = NULL;
PDT *DT = NULL;
std::vector<Pair_P_I> pos_idx_pair;
std::vector<Bird_wo_pos> birds;

void ini_rand(double _Lx, double _Ly, double _eta, unsigned long long seed, int N) {
  Lx = _Lx;
  Ly = _Ly;
  eta = _eta;
  nStep = N;
  myran = new Ran(seed);
  unsigned int nBird = int (rho0 * Lx * Ly);

	// initialize birds randomly
	Bird_wo_pos::ini_rand(birds, pos_idx_pair, nBird, Lx, Ly, v0, *myran);

	// initialize triangulation
  PDT::Iso_rectangle domain(0, 0, Lx, Ly);
  DT = new PDT(domain);

  std::cout << "number of birds: " << nBird << "\n";
  std::cout << "system size: " << Lx << ", " << Ly << "\n";
  std::cout << "eta = " << eta << "\n";
}

void run() {
	align(*DT, pos_idx_pair, birds);
	move_with_vectorial_noise(eta, *myran, pos_idx_pair, birds);
}

void run(int n) {
  for (int i = 0; i < n; i++){
    align(*DT, pos_idx_pair, birds);
    move_with_vectorial_noise(eta, *myran, pos_idx_pair, birds);
  }
}

void run(int n, double &phi) {
	run(n);
	phi = cal_order_para(birds);
}

void run_test(int n) {
   for (int i = 0; i < n; i++){
    align(*DT, pos_idx_pair, birds);
    move_with_vectorial_noise(eta, *myran, pos_idx_pair, birds);
		double phi = cal_order_para(birds);
    std::cout << i << "\t" << phi << "\n";
  } 
}
