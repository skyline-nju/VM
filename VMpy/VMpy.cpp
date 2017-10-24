#include "VMpy.h"
#include <iostream>
#ifdef _MSC_VER
#include "../VM/run.h"
#else
#include "run.h"
#endif

Ran *myran = NULL;
Node *bird = NULL;
Grid *cell = NULL;
double *disorder = NULL;
double eta;
double eps;

void set_random_seed(int seed) {
	myran = new Ran(seed);
}

void set_v0(double _v0) {
	Node::v0 = _v0;
}

void set_eta(double _eta) {
	eta = _eta;
}

void set_eps(double _eps) {
	eps = _eps;
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

void ini(double * x, double * y, double * vx, double * vy, int nBird,
				 int seed, double v0, double _eta, double _eps, double Lx, double Ly) {
	myran = new Ran(seed);
	ini_rand_torques(&disorder, Lx * Ly, _eps, seed);
	Node::v0 = v0;
	eta = _eta;
	eps = _eps;
	Node::Lx = Lx;
	Node::Ly = Ly;
	Node::N = nBird;
	Node::rho_0 = nBird / (Node::Lx * Node::Ly);

	cell = Grid::ini(Node::Lx, Node::Ly);
	bird = new Node[Node::N];
	pack(x, y, vx, vy, nBird, bird);
	Grid::link_nodes(cell, bird);
}

void run(int nstep) {
	run_raw(bird, cell, myran, nstep, eta, NULL);
}

void get_snap(double * x, double * y, double * vx, double * vy, int nBird) {
	unpack(x, y, vx, vy, nBird, bird);
}

void get_coarse_grained_snap(int *num, double * svx, double * svy,
														 int ncells, double l) {
	for (int i = 0; i < ncells; i++) {
		num[i] = 0;
		svx[i] = 0;
		svy[i] = 0;
	}
	int ncols = int(Node::Lx / l);
	for (int i = 0; i < Node::N; i++) {
		int ic = int(bird[i].x / l) + ncols * int(bird[i].y / l);
		svx[ic] += bird[i].vx;
		svy[ic] += bird[i].vy;
		num[ic] += 1;
	}
}

