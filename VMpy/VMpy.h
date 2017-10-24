#ifndef VMPY_H
#define VMPY_H

void set_random_seed(int seed);
void set_v0(double _v0);
void set_eta(double _eta);
void set_eps(double _eps);
void setLx(double _Lx);
void setLy(double _Ly);

void ini(double *x, double *y, double *vx, double *vy, int nBird,
				 int seed, double v0, double _eta, double _eps, double Lx, double Ly);

void run(int nstep);

void get_snap(double *x, double *y, double *vx, double *vy, int nBird);

void get_coarse_grained_snap(int *num, double *svx, double *svy,
														 int ncells, double l);
#endif
