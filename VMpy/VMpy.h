#ifndef VMPY_H
#define VMPY_H

void set_random_seed(int seed);
void set_v0(double _v0);
void set_eta(double _eta);
void setLx(double _Lx);
void setLy(double _Ly);

void ini(double *x, double *y, double *vx, double *vy, int nBird,
				 int seed, double v0, double _eta, double Lx, double Ly);

void run(int nstep);

void get_snap(double *x, double *y, double *vx, double *vy, int nBird);

void get_coarse_grained_snap(int *num, double *svx, double *svy,
	                           int ncells, double l);
#endif
