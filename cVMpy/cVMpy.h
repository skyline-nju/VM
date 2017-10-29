#ifndef CVMPY_H
#define CVMPY_H

void ini(double *x, double *y, double *theta, int nBird, double _eta,
				 double Lx, double Ly, int seed, double dt);

void ini_rand(int nBird, double _eta, double Lx, double Ly, int seed, double dt);

void update_coor();

void run(int nstep);

void get_snap(double *x, double *y, double *theta, int nBird);

void set_h(double _h);

void set_period(int _period);

double get_h();

#endif
