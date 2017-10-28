#ifndef XYMODEL_H
#define XYMODEL_H

const double PI = 3.14159265358979;

void ini(int _L, double _dt, double _eta, int seed,
					double *_phi, int n_phi);

void ini(int _L, double _dt, double _eta, int seed,
					double *_phi, int n_phi, double _h, double *_psi);

void one_step();

void run(int nstep, double *phi_out, int size_phi);
#endif
