#include "XYmodel.h"
#ifdef _MSC_VER
#include "../VM/rand.h"
#else
#include "rand.h"
#endif
#include <cmath>
#include <cstring>

int L;
int n;
double dt;
double noise_strength;
double h;						  // Strength of the external fields
double *phi;
double *phi_next;
double *psi=NULL;
int *neighbor;
Ran *myran;

void set_neighbor() {
	neighbor = new int[n * 8];
	int dcol[] = { -1, 0, 1, -1, 1, -1, 0, 1 };
	int drow[] = { -1, -1, -1, 0, 0, 1, 1, 1 };
	for (int i = 0; i < n; i++) {
		int col_i = i % L;
		int row_i = i / L;
		for (int j = 0; j < 8; j++) {
			int col_j = col_i + dcol[j];
			int row_j = row_i + drow[j];
			if (col_j < 0) {
				col_j += L;
			} else if (col_j >= L) {
				col_j -= L;
			}
			if (row_j < 0) {
				row_j += L;
			} else if (row_j >= L) {
				row_j -= L;
			}
			neighbor[i * 8 + j] = col_j + row_j * L;
		}
	}
}

void ini(int _L, double _dt, double _eta, int seed, 
	                double *_phi, int n_phi) {
	L = _L;
	dt = _dt;
	noise_strength = _eta * 2 * PI;
	n = n_phi;
	phi = new double[n];
	phi_next = new double[n];
	myran = new Ran(seed);
	for (int i = 0; i < n; i++) {
		phi[i] = _phi[i];
	}
	set_neighbor();
}

void ini(int _L, double _dt, double _eta, int seed,
					double * _phi, int n_phi, double _h, double *_psi) {
	ini(_L, _dt, _eta, seed, _phi, n_phi);
	h = _h;
	psi = new double[n_phi];
	memcpy(psi, _psi, sizeof(double) * n_phi);
}

void one_step() {
	for (int i = 0; i < n; i++) {
		double noise = noise_strength * (myran->doub() - 0.5);
		int j0 = i * 8;
		for (int j = 0; j < 8; j++) {
			noise += sin(phi[i] - phi[neighbor[j0 + j]]);
		}
		if (psi) {
			noise += sin(phi[i] - psi[i]) * h;
		}
		phi_next[i] = phi[i] - dt * noise;
	}

	double *tmp = phi;
	phi = phi_next;
	phi_next = tmp;
}

void run(int nstep, double *phi_out, int size_phi) {
	for (int i = 0; i < nstep; i++) {
		one_step();
	}
	memcpy(phi_out, phi, sizeof(double) * size_phi);
}
