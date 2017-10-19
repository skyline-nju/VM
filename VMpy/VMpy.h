#ifndef VMPY_H
#define VMPY_H

void set_random_seed(int seed);
void set_v0(double _v0);
void set_eta(double _eta);
void setLx(double _Lx);
void setLy(double _Ly);

void run(int nstep, double *x, double *y, double *vx, double *vy, int nBird);

void coarse_grain(double l, double *theta, int ncells,
									double *x, double *y, double *vx, double *vy, int nBird);
#endif
