#pragma once
#include "config.h"
#include <cstdlib>
#include <vector>

struct Node {
	Node() : x(0.), y(0.), vx(0.), vy(0.), align_x(0.), align_y(0.), replusion_x(0.), replusion_y(0.), n_neighbor(0), next(NULL) {}
	double rr(Node* node) {
		double dx = node->x - x;
		double dy = node->y - y;
		return dx * dx + dy * dy;
	}
	double rr(Node* node, double a, double b) {
		double dx = node->x - x + a;
		double dy = node->y - y + b;
		return dx * dx + dy * dy;
	}

	double cal_dis2(const Node* p2, double& dx, double& dy) {
		dx = p2->x - x;
		dy = p2->y - y;
		return dx * dx + dy * dy;
	}

	double cal_dis2(const Node* p2, double& dx, double& dy, double a, double b) {
		dx = p2->x - x + a;
		dy = p2->y - y + b;
		return dx * dx + dy * dy;
	}

	void align_repel(Node* p2, double dr_hat_x, double dr_hat_y);

	void align_repel(Node* p2, double dr_hat_x, double dr_hat_y, double beta_h, double beta_v);

	void collide(Node* p2);

	void collide(Node* p2, double a, double b);

	double x;
	double y;
	double vx;
	double vy;
	double align_x;
	double align_y;
	double replusion_x;
	double replusion_y;
	int n_neighbor;
	Node* next;

#ifdef MULT_REP
	static double beta_h;
	static double beta_v;
#endif
};

struct Coordinate {
	float x;
	float y;
	float theta;
};


//void alignment(Grid* grid);