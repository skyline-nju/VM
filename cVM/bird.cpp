#include "bird.h"

int Bird::N;
double Bird::rho0;
double Bird::v0;
double Bird::dt;
double Bird::Lx;
double Bird::Ly;

void Bird::interact(Bird *bird) {
	if (rr(bird) < 1) {
		double s = sin(bird->theta - theta);
		omega_align += s;
		bird->omega_align -= s;
		tot_neighbor++;
		bird->tot_neighbor++;
	}
}

void Bird::interact(Bird *bird, double a, double b) {
	if (rr(bird, a, b) < 1) {
		double s = sin(bird->theta - theta);
		omega_align += s;
		bird->omega_align -= s;
		tot_neighbor++;
		bird->tot_neighbor++;
	}
}

void Bird::move(double omega_ext) {
	theta += (omega_align + omega_ext) * dt;
	omega_align = 0;
	x += v0 * dt * cos(theta);
	y += v0 * dt * sin(theta);

	if (x >= Lx)
		x -= Lx;
	else if (x < 0)
		x += Lx;
	if (y >= Ly)
		y -= Ly;
	else if (y < 0)
		y += Ly;
}

void Bird::move(double omega_ext, double h) {
	if (tot_neighbor > 0)
		omega_align /= tot_neighbor;
	theta += (omega_align + omega_ext - h* sin(theta)) * dt;
	omega_align = 0;
	tot_neighbor = 0;
	x += v0 * dt * cos(theta);
	y += v0 * dt * sin(theta);

	if (x >= Lx)
		x -= Lx;
	else if (x < 0)
		x += Lx;
	if (y >= Ly)
		y -= Ly;
	else if (y < 0)
		y += Ly;
}
