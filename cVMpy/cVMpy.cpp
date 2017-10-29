#include "cVMpy.h"
#ifdef _MSC_VER
#include "../cVM/rand.h"
#include "../cVM/bird.h"
#include "../cVM/cell_list.h"
#else
#include "rand.h"
#include "bird.h"
#include "cell_list.h"
#endif
#include <iostream>

const double PI = 3.14159265358979;
Bird *bird = NULL;
Cell<Bird> *cell = NULL;
Ran *myran = NULL;
double eta;
double h = 0;
int period = 0;

void set_h(double _h) {
	h = _h;
}

double get_h() {
	return h;
}
void set_period(int _period) {
	period = _period;
}


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
		//bird[i].move(noise);
		if (h != 0)
			bird[i].move(noise, h);
		else
			bird[i].move(noise);
	}
}

void run(int nstep) {
	for (int i = 1; i <= nstep; i++) {
		Cell<Bird>::all_pairs(cell);
		update_coor();
		Cell<Bird>::refresh(cell, bird);
		if (period > 0 && i % period == 0) {
			h *= -1;
		}
	}
}

void get_snap(double *x, double *y, double *theta, int nBird) {
	for (int i = 0; i < nBird; i++) {
		x[i] = bird[i].x;
		y[i] = bird[i].y;
		theta[i] = bird[i].theta;
	}
}

