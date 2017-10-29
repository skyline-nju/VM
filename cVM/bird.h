#ifndef BIRD_H
#define BIRD_H
#include <cmath>
#include <cstdio>

class Bird{
public:
	Bird() { next = NULL; tot_neighbor = 0; }
	double rr(Bird *node);
	double rr(Bird *node, double a, double b);
	void interact(Bird *bird);
	void interact(Bird *bird, double a, double b);
	void move(double omega_ext);
	void move(double omega_ext, double h);

	double x;
	double y;
	//double ux;
	//double uy;
	double theta;	  
	double omega_align;  // angular velocity due to alignment
	int cell_idx;
	int tot_neighbor;
	Bird *next;

	static int N;
	static double rho0;
	static double v0;
	static double dt;
	static double Lx;
	static double Ly;
};

inline double Bird::rr(Bird *bird) {
	double dx = bird->x - x;
	double dy = bird->y - y;
	return dx*dx + dy*dy;
}

inline double Bird::rr(Bird *bird, double a, double b) {
	double dx = bird->x - x + a;
	double dy = bird->y - y + b;
	return dx*dx + dy*dy;
}
#endif
