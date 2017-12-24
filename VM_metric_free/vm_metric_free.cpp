#include <fstream>
#include <ctime>
#include <chrono>
#include "bird.h"

// calculate the instant order parameter
double cal_order_para(std::vector<Bird_wo_pos> &bird) {
	double svx = 0;
	double svy = 0;
	for (auto it = bird.cbegin(); it != bird.cend(); ++it) {
		svx += (*it).vx;
		svy += (*it).vy;
	}
	double phi = sqrt(svx * svx + svy * svy) / bird.size();
	return phi;
}

/******************************** MAIN *************************************/
int main(int argc, char* argv[]) {
	// set parameters
	double Lx = atof(argv[1]);
	double Ly = Lx;
	double eta = atof(argv[2]);
	Ran myran(atoi(argv[3]));
	int nstep = atoi(argv[4]);
	double rho0 = 1;
	double v0 = 0.5;
	
	unsigned int nBird = unsigned int(rho0 * Lx * Ly);

	// initialize birds randomly
	std::vector<Pair_P_I> pos_idx_pair;
	std::vector<Bird_wo_pos> birds;
	Bird_wo_pos::ini_rand(birds, pos_idx_pair, nBird, Lx, Ly, v0, myran);

	// initialize triangulation
  PDT::Iso_rectangle domain(0, 0, Lx, Ly);
  PDT DT(domain);

	std::cout << "number of birds: " << nBird << "\n";
	std::cout << "system size: " << Lx << ", " << Ly << "\n";

	for (int step = 0; step < nstep; step++) {
		if (step % 100 == 0)
			std::cout << step << "\t" << cal_order_para(birds) << "\n";
		align(DT, pos_idx_pair, birds);
		move_with_vectorial_noise(eta, myran, pos_idx_pair, birds);
	}
}
