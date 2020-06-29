#include "particle.h"
#include <cmath>

#ifdef MULT_REP
	double Node::beta_h;
	double Node::beta_v;
#endif
void Node::align_repel(Node* p2, double dr_hat_x, double dr_hat_y) {
	align_x += p2->vx;
	align_y += p2->vy;
	p2->align_x += vx;
	p2->align_y += vy;
#ifdef ISO_REP
	replusion_x -= dr_hat_x;
	replusion_y -= dr_hat_y;
	p2->replusion_x += dr_hat_x;
	p2->replusion_y += dr_hat_y;
#else
#ifdef ANI_REP_0
	double rep1 = vx * dr_hat_x + vy * dr_hat_y;
	double rep2 = -p2->vx * dr_hat_x - p2->vy * dr_hat_y;
#elif defined ANI_REP_90
	double rep1 = vx * dr_hat_y - vy * dr_hat_x;
	double rep2 = -p2->vx * dr_hat_y + p2->vy * dr_hat_x;
#elif defined ANI_REP_0_90
	double rep1 = (vx * dr_hat_x + vy * dr_hat_y) * (vx * dr_hat_y - vy * dr_hat_x);
	double rep2 = (-p2->vx * dr_hat_x - p2->vy * dr_hat_y) * (-p2->vx * dr_hat_y + p2->vy * dr_hat_x);
#endif
	rep1 = rep1 * rep1;
	rep2 = rep2 * rep2;
	replusion_x -= rep1 * dr_hat_x;
	replusion_y -= rep1 * dr_hat_y;
	p2->replusion_x += rep2 * dr_hat_x;
	p2->replusion_y += rep2 * dr_hat_y;
	//double rep_m = 0.5 * (rep1 + rep2);
	//replusion_x -= rep_m * dr_hat_x;
	//replusion_y -= rep_m * dr_hat_y;
	//p2->replusion_x += rep_m * dr_hat_x;
	//p2->replusion_y += rep_m * dr_hat_y;
#endif
	n_neighbor++;
	p2->n_neighbor++;
}

void Node::align_repel(Node* p2, double dr_hat_x, double dr_hat_y, double beta_h, double beta_v) {
	align_x += p2->vx;
	align_y += p2->vy;
	p2->align_x += vx;
	p2->align_y += vy;

	double rep1_v = vx * dr_hat_x + vy * dr_hat_y;
	double rep2_v = -p2->vx * dr_hat_x - p2->vy * dr_hat_y;

	double rep1_h = vx * dr_hat_y - vy * dr_hat_x;
	double rep2_h = -p2->vx * dr_hat_y + p2->vy * dr_hat_x;

	double rep1 = beta_v * rep1_v * rep1_v + beta_h * rep1_h * rep1_h;
	double rep2 = beta_v * rep2_v * rep2_v + beta_h * rep2_h * rep2_h;
	replusion_x -= rep1 * dr_hat_x;
	replusion_y -= rep1 * dr_hat_y;
	p2->replusion_x += rep2 * dr_hat_x;
	p2->replusion_y += rep2 * dr_hat_y;
	//double rep_m = 0.5 * (rep1 + rep2);
	//replusion_x -= rep_m * dr_hat_x;
	//replusion_y -= rep_m * dr_hat_y;
	//p2->replusion_x += rep_m * dr_hat_x;
	//p2->replusion_y += rep_m * dr_hat_y;

	n_neighbor++;
	p2->n_neighbor++;
}

void Node::collide(Node* p2) {
	double dx, dy;
	double rr = cal_dis2(p2, dx, dy);
	if (rr < 1.) {
		double inv_r = 1. / sqrt(rr);
		double dr_hat_x = dx * inv_r;
		double dr_hat_y = dy * inv_r;
#ifndef MULT_REP
		align_repel(p2, dr_hat_x, dr_hat_y);
#else
		align_repel(p2, dr_hat_x, dr_hat_y, beta_h, beta_v);
#endif
	}
}

void Node::collide(Node* p2, double a, double b) {
	double dx, dy;
	double rr = cal_dis2(p2, dx, dy, a, b);
	if (rr < 1.) {
		double inv_r = 1. / sqrt(rr);
		double dr_hat_x = dx * inv_r;
		double dr_hat_y = dy * inv_r;
#ifndef MULT_REP
		align_repel(p2, dr_hat_x, dr_hat_y);
#else
		align_repel(p2, dr_hat_x, dr_hat_y, beta_h, beta_v);
#endif
	}
}