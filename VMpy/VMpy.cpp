#include "VMpy.h"
#include <iostream>
#include "../VM/node.h"
#include "../VM/grid.h"
#include "../VM/run.h"

// initialize random number generator
Ran *myran = NULL;
int my_seed;
double eta;

void set_random_seed(int seed) {
	myran = new Ran(seed);
	my_seed = seed;
}

void set_v0(double _v0) {
	Node::v0 = _v0;
}

void set_eta(double _eta) {
	eta = _eta;
}

void setLx(double _Lx) {
	Node::Lx = _Lx;
}

void setLy(double _Ly) {
	Node::Ly = _Ly;
}

void pack(const double *x, const double *y, const double *vx, const double *vy,
					int nBird, Node * bird) {
	for (int i = 0; i < nBird; i++) {
		bird[i].x = x[i];
		bird[i].y = y[i];
		bird[i].vx = bird[i].vx0 = vx[i];
		bird[i].vy = bird[i].vy0 = vy[i];
		bird[i].next = NULL;
 }
}

void unpack(double *x, double *y, double *vx, double *vy,
						int nBird, const Node * bird) {
	for (int i = 0; i < nBird; i++) {
		x[i] = bird[i].x;
		y[i] = bird[i].y;
		vx[i] = bird[i].vx;
		vy[i] = bird[i].vy;
	}
}

void run(int nstep, double *x, double *y, double *vx, double *vy, int nBird) {
	Node::N = nBird;
	Node::rho_0 = nBird / (Node::Lx * Node::Ly);

	// initialize cell list; calculate total number of cells.
	Grid *cell = Grid::ini(Node::Lx, Node::Ly);
	
	// initialize location and velocity of birds
	Node *bird = new Node[nBird];
	pack(x, y, vx, vy, nBird, bird);

	// link birds to cell list
	Grid::link_nodes(cell, bird);

	// run nstep time steps
	run_raw(bird, cell, myran, nstep, eta, NULL);

	// unpack struct to arrays
	unpack(x, y, vx, vy, nBird, bird);

	delete[] bird;
	delete[] cell;
}

void coarse_grain(double l, double *theta, int ncells,
									double *x, double *y, double *vx, double *vy, int nBird) {
	double *svx = new double[ncells];
	double *svy = new double[ncells];
	for (int i = 0; i < ncells; i++) {
		theta[i] = 0;
		svx[i] = 0;
		svy[i] = 0;
	}
	int ncols = int(Node::Lx / l);
	int nrows = int(Node::Ly / l);
	for (int i = 0; i < nBird; i++) {
		int ic = int(x[i] / l) + ncols * int(y[i] / l);
		svx[ic] += vx[i];
		svy[ic] += vy[i];
	}
	for (int i = 0; i < ncells; i++) {
		theta[i] = atan2(svy[i], svx[i]);
	}
}

