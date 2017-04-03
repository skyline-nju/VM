#include "run.h"
#include "iodata.h"
using namespace std;

void ini_rand_torques(double **disorder, int n, double epsilon, Ran *myran)
{
	double d = 1.0 / (n - 1);
	double *torque = new double[n];
	for (int i = 0; i < n; i++)
		torque[i] = (-0.5 + i *d) * epsilon * 2.0 * PI;
	shuffle(torque, n, myran);
	*disorder = torque;
}

void update_coor(Node *bird, Ran *myran, double eta)
{
	static double eta2PI = eta * 2 * PI;
	double *noise = new double[Node::N];
	for (int i = 0; i < Node::N; i++)
	{
		noise[i] = (myran->doub() - 0.5) * eta2PI;
	}
	for (int i = 0; i < Node::N; i++)
	{
		bird[i].move(noise[i]);
	}
	delete[] noise;
	noise = nullptr;
}

void update_coor(Node *bird, Ran* myran, double eta, double epsilon, const double *disorder)
{
	static double eta2PI = eta * 2 * PI;

	//calculating noise
	double *noise = new double[Node::N];

	for (int i = 0; i < Node::N; i++)
	{
		noise[i] = (myran->doub() - 0.5) * eta2PI + disorder[bird[i].cell_idx];
	}

	//updating coordination
	for (int i = 0; i < Node::N; i++)
	{
		bird[i].move(noise[i]);
	}
	delete[] noise;
	noise = nullptr;
}

void run(Node *bird, Grid *cell, Ran*myran, int nStep, double eta)
{
	for (int i = 1; i <= nStep; i++)
	{
		Grid::all_pairs(cell);
		update_coor(bird, myran, eta);
		output(bird, i);
	}
}

void run(Node *bird, Grid *cell, Ran*myran, int nStep, double eta, double epsilon, const double *disorder)
{
	for (int i = 1; i <= nStep; i++)
	{
		Grid::all_pairs(cell);
		update_coor(bird, myran, eta, epsilon, disorder);
		output(bird, i);
	}
}
