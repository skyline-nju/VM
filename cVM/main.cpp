#include "rand.h"
#include "bird.h"
#include "cell_list.h"
#include <iostream>

const double PI = 3.14159265358979;
Bird *bird = NULL;
Cell<Bird> *cell = NULL;
Ran *myran = NULL;
double eta;

void ini(double * x, double * y, double * theta, int nBird, double _eta,
				 double Lx, double Ly, int seed, double dt) {
	// initailize random number generator
	myran = new Ran(seed);

	// initial location and orientation of birds
	bird = new Bird[nBird];
	for (int i = 0; i < nBird; i++) {
		bird[i].x = x[i];
		bird[i].y = y[i];
		bird[i].theta = theta[i];
		bird[i].omega_align = 0;
		bird[i].next = NULL;
	}
	Bird::N = nBird;
	Bird::rho0 = nBird / (Lx * Ly);
	Bird::v0 = 0.5;
	Bird::Lx = Lx;
	Bird::Ly = Ly;
	Bird::dt = dt;
	eta = _eta;

	// initailize cell list
	Cell<Bird>::l0 = 1;
	cell = Cell<Bird>::ini(bird, Lx, Ly);
}

void ini_rand(int nBird, double _eta, double Lx, double Ly, int seed, double dt) {
	// initailize random number generator
	myran = new Ran(seed);

	// initial location and orientation of birds
	bird = new Bird[nBird];
	for (int i = 0; i < nBird; i++) {
		bird[i].x = myran->doub() * Lx;
		bird[i].y = myran->doub() * Ly;
		bird[i].theta = myran->doub() * PI * 2;
		bird[i].omega_align = 0;
		bird[i].next = NULL;
	}
	Bird::N = nBird;
	Bird::rho0 = nBird / (Lx * Ly);
	Bird::v0 = 0.5;
	Bird::Lx = Lx;
	Bird::Ly = Ly;
	Bird::dt = dt;
	eta = _eta;

	// initailize cell list
	Cell<Bird>::l0 = 1;
	cell = Cell<Bird>::ini(bird, Lx, Ly);

}

void update_coor() {
	static double eta2PI = eta * 2 * PI;
	for (int i = 0; i < Bird::N; i++) {
		double noise = (myran->doub() - 0.5) * eta2PI;
		bird[i].move(noise);
	}
}

void run(int nstep) {
	for (int i = 1; i <= nstep; i++) {
		Cell<Bird>::all_pairs(cell);
		update_coor();
		Cell<Bird>::refresh(cell, bird);
	}
}

void get_snap(double *x, double *y, double *theta, int nBird) {
	for (int i = 0; i < nBird; i++) {
		x[i] = bird[i].x;
		y[i] = bird[i].y;
		theta[i] = bird[i].theta;
	}
}

int main() {
	ini_rand(100, 0.1, 100, 100, 123, 0.1);
	run(100);
	std::cout << "hello, world" << std::endl;
}