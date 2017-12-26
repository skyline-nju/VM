#include <fstream>
#include <ctime>
#include <chrono>
#include "metric_free_vm.h"
#include "comn.h"

#ifdef _MSC_VER
std::string delimiter("\\");
#else
std::string delimiter("/");
#endif

int main(int argc, char* argv[]) {
	// set parameters
	double Lx = atof(argv[1]);
	double Ly = Lx;
	double eta = atof(argv[2]);
	unsigned long long seed = atoi(argv[3]);
	int nstep = atoi(argv[4]);
	double eps = 0;

	ini_rand(Lx, Ly, eta, seed, nstep);

	// output setting
	mkdir("phi");
	char filename[100];
	snprintf(filename, 100, "phi%s%g_%g_%g_%llu.dat", delimiter.c_str(), Lx, eta, eps, seed);
	std::ofstream fout_phi(filename);
	//mkdir("log");
	//snprintf(filename, 100, "log%s%g_%g_%g_%llu.dat", delimiter.c_str(), Lx, eta, 0, seed);
	int dt_phi = 100;
	int n_period = nstep / dt_phi;

	auto t_start = std::chrono::system_clock::now();
	for (int i = 0; i < n_period; i++) {
		double phi;
		run(dt_phi, phi);
		auto t_now = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_time = t_now - t_start;
		fout_phi << (i + 1) * dt_phi << "\t" << phi << "\t" << elapsed_time.count() << "\n";
	}
}
