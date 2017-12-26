#include "bird.h"
#include "comn.h"

double Bird_wo_pos::v0;
double Bird_wo_pos::Lx;
double Bird_wo_pos::Ly;

void Bird_wo_pos::move(Point &p) const {
	double x = p[0] + v0 * vx;
	double y = p[1] + v0 * vy;
	if (x >= Lx) {
		x -= Lx;
	} else if (x < 0) {
		x += Lx;
	}
	if (y >= Ly) {
		y -= Ly;
	} else if (y < 0) {
		y += Ly;
	}
	p += Vec_2(x-p[0], y-p[1]);  
}

void Bird_wo_pos::vectorial_noise(double eta, Ran & myran) {
	double noise_x, noise_y;
	myran.circle_point_picking(noise_x, noise_y);
	double tmp = eta * n_neighbor;
	vx_next += noise_x * tmp;
	vy_next += noise_y * tmp;
	tmp = 1 / sqrt(vx_next * vx_next + vy_next * vy_next);
	vx_next *= tmp;
	vy_next *= tmp;
	vx = vx_next;
	vy = vy_next;
	n_neighbor = 1;
}

void Bird_wo_pos::scalar_noise(double eta, Ran & myran) {
	double tmp = std::sqrt(vx_next * vx_next + vy_next * vy_next);
	double c1 = vx / tmp;
	double s1 = vy / tmp;
	double noise = myran.doub() * eta * 2 * PI;
	double c2 = cos(noise);
	double s2 = sin(noise);
	vx = vx_next = c1 * c2 - s1 * s2;
	vy = vy_next = c1 * s2 + c2 * s1;
	n_neighbor = 1;
}



