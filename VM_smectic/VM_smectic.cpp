/*********************************************************************
Vicsek model with smectic order
2020/6/14
*********************************************************************/
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <cstring>
#include <cstdio>
#include "rand.h"
#include "config.h"
#include "particle.h"
#include "exporter.h"
using namespace std;
ofstream fout;

double L, rho, l0 = 1.0, R = 1.0, RR = R * R, v0 = 0.3;
double sigma, sigma2PI, beta_h, beta_v;
int N, m, mm, stp;
int seed;
Ran* myran = NULL;

struct Grid
{
	Grid() { head = NULL; }
	void addNode(Node* node) {
		node->next = head;
		head = node;
	}
	void align();
	void align(Grid*);
	void align(Grid*, double, double);
	Node* head;
};

void Grid::align()
{
	Node* node1 = head;
	Node* node2;
	while (node1->next)
	{
		node2 = node1->next;
		do {
			node1->collide(node2);
			node2 = node2->next;
		} while (node2);
		node1 = node1->next;
	}
}

void Grid::align(Grid* grid)
{
	if (grid->head)
	{
		Node* node1 = head;
		Node* node2;
		do
		{
			node2 = grid->head;
			do {
				node1->collide(node2);
				node2 = node2->next;
			} while (node2);
			node1 = node1->next;
		} while (node1);
	}
}

void Grid::align(Grid* grid, double a, double b)
{
	if (grid->head)
	{
		Node* node1 = head;
		Node* node2;
		do
		{
			node2 = grid->head;
			do {
				node1->collide(node2, a, b);
				node2 = node2->next;
			} while (node2);
			node1 = node1->next;
		} while (node1);
	}
}
void alignment(Grid* grid)
{
	int i, j;
	Grid* p = grid;
	for (j = 0; j <= m - 2; j++)
	{
		if (p->head)
		{
			p->align();
			p->align(p + 1);
			p->align(p + m + m - 1, -L, 0);
			p->align(p + m);
			p->align(p + m + 1);
		}
		p++;
		for (i = 1; i <= m - 2; i++)
		{
			if (p->head)
			{
				p->align();
				p->align(p + 1);
				p->align(p + m - 1);
				p->align(p + m);
				p->align(p + m + 1);
			}
			p++;
		}
		if (p->head)
		{
			p->align();
			p->align(p - m + 1, L, 0);
			p->align(p + 1, L, 0);
			p->align(p + m);
			p->align(p + m - 1);
		}
		p++;
	}
	if (p->head)
	{
		p->align();
		p->align(p + 1, 0, 0);
		p->align(grid, 0, L);
		p->align(grid + 1, 0, L);
		p->align(grid + m - 1, -L, L);
	}
	p++;
	for (i = 1; i <= m - 2; i++)
	{
		if (p->head)
		{
			p->align();
			p->align(p + 1);
			p->align(grid + i - 1, 0, L);
			p->align(grid + i, 0, L);
			p->align(grid + i + 1, 0, L);
		}
		p++;
	}
	if (p->head)
	{
		p->align();
		p->align(p - m + 1, L, 0);
		p->align(grid + m - 2, 0, L);
		p->align(grid + m - 1, 0, L);
		p->align(grid, L, L);
	}
}

void ini_rand(Node* bird, Grid* grid) {
	for (int i = 0; i < N; i++) {
		double x = myran->doub() * L;
		double y = myran->doub() * L;
		double theta = myran->doub() * PI * 2.0;
		bird[i].x = x;
		bird[i].y = y;
		bird[i].vx = cos(theta);
		bird[i].vy = sin(theta);
		bird[i].align_x = bird[i].align_y = 0.;
		bird[i].replusion_x = bird[i].replusion_y = 0.;
		bird[i].n_neighbor = 0;
		bird[i].next = NULL;
		int ic = int(x) + int(y) * m;
		grid[ic].addNode(bird + i);
	}
}

void move_forward(Node* bird, Grid* grid) {
	for (int ic = 0; ic < mm; ic++) {
		grid[ic].head = NULL;
	}

	for (int i = 0; i < N; i++) {
		double x = bird[i].x;
		double y = bird[i].y;
		int tot_neighbor = bird[i].n_neighbor + 1;
		double c1, s1;
		if (tot_neighbor > 1) {
			c1 = (bird[i].vx + bird[i].align_x) / tot_neighbor + bird[i].replusion_x / bird[i].n_neighbor;
			s1 = (bird[i].vy + bird[i].align_y) / tot_neighbor + bird[i].replusion_y / bird[i].n_neighbor;
		} else {
			c1 = bird[i].vx + bird[i].align_x;
			s1 = bird[i].vy + bird[i].align_y;
		}

		double tmp = sqrt(c1 * c1 + s1 * s1);
		c1 /= tmp;
		s1 /= tmp;
		double noise = sigma2PI * (myran->doub() - 0.5);
		double c2 = cos(noise);
		double s2 = sin(noise);
		double vx = c1 * c2 - s1 * s2;
		double vy = c1 * s2 + c2 * s1;

		x += v0 * vx;
		if (x >= L) x -= L;
		else if (x < 0) x += L;
		y += v0 * vy;
		if (y >= L) y -= L;
		else if (y < 0) y += L;

		bird[i].x = x;
		bird[i].y = y;
		bird[i].vx = vx;
		bird[i].vy = vy;
		bird[i].align_x = 0.;
		bird[i].align_y = 0.;
		bird[i].replusion_x = 0.;
		bird[i].replusion_y = 0.;
		bird[i].n_neighbor = 0;
	}

	for (int i = 0; i < N; i++) {
		int ic = int(bird[i].x) + int(bird[i].y) * m;
		//if (ic < 0 || ic >= mm) {
		//	std::cout << bird[i].x << ", " << bird[i].y << std::endl;
		//	std::cout << "ic = " << ic << std::endl;
		//}
		Node* p = &bird[i];
		grid[ic].addNode(p);
	}
}

void order(Node* bird) {
	double svx = 0, svy = 0, phi, theta;
	for (int i = 0; i < N; i++) {
		svx += bird[i].vx;
		svy += bird[i].vy;
	}
	phi = sqrt(svx * svx + svy * svy) / N;
	theta = atan2(svy, svx);
	fout << fixed << setw(16) << setprecision(10) << phi << "\t" << theta << endl;
}

int main(int argc, char* argv[]) {
	sigma = atof(argv[1]);
	beta_h = atof(argv[2]);
	beta_v = atof(argv[3]);
	rho = atof(argv[4]);
	L = atof(argv[5]);
	stp = atoi(argv[6]);
	seed = atoi(argv[7]);

	myran = new Ran(seed);
	sigma2PI = sigma * 2.0 * PI;

	Node::beta_h = beta_h;
	Node::beta_v = beta_v;

	char file[100];
	snprintf(file, 100, "s%.2f_b%.2f_%.2f_r%.2f_L%g_v%g.extxyz", sigma, beta_h, beta_v, rho, L, v0);
	exporter::XyzExporter_2 xy_exp(file, 0, stp, 1000, Vec_2<double>(L, L));
	snprintf(file, 100, "s%.2f_b%.2f_%.2f_r%.2f_L%g_v%g.dat", sigma, beta_h, beta_v, rho, L, v0);
	fout.open(file);

	N = int(L * L * rho);
	m = int(L / l0);
	mm = m * m;

	cout << "L=" << L << endl;
	cout << "m=" << m << endl;
	cout << "sigma=" << sigma << endl;
	cout << "beta=" << beta_h << "\t" << beta_v << endl;
	cout << "rho=" << rho << endl;
	cout << "Nstep=" << stp << endl;
	cout << "seed=" << seed << endl;
	cout << endl;
	Node* bird = new Node[N];
	Grid* grid = new Grid[mm];

	ini_rand(bird, grid);

	std::cout << "initialized!" << std::endl;
	//clock_t t_beg = clock();

	xy_exp.dump_pos_ori(0, bird, N);
	for (int j = 1; j <= stp; j++) {
		//std::cout << "j = " << j << std::endl;
		alignment(grid);
		move_forward(bird, grid);
		xy_exp.dump_pos_ori(j, bird, N);
		if (j % 100 == 0) {
			order(bird);
		}
	}
	fout.close();
	return 0;
}