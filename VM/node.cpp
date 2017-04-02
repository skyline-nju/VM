#include "node.h"

using namespace std;

double Node::Lx;
double Node::Ly;
double Node::v0 = 0.5;
double Node::rho_0;
int Node::N;

void Node::align(Node *node)
{
	if (rr(node) < 1)
		addV(node);
}

void Node::align(Node *node, double a, double b)
{
	if (rr(node, a, b) < 1)
		addV(node);
}

void Node::move(double noise)
{
	double tmp = sqrt(vx*vx + vy*vy);
	double c1 = vx / tmp;
	double s1 = vy / tmp;
	double c2 = cos(noise);
	double s2 = sin(noise);
	vx = vx0 = c1 * c2 - s1 * s2;
	vy = vy0 = c1 * s2 + c2 * s1;

	x += v0*vx;
	if (x >= Lx)
		x -= Lx;
	else if (x < 0)
		x += Lx;
	y += v0*vy;
	if (y >= Ly)
		y -= Ly;
	else if (y < 0)
		y += Ly;

}

Node * Node::ini_rand(Ran * myran)
{
	N = int(rho_0 * Lx * Ly);
	Node *bird = new Node[N];
	for (int i = 0; i < N; i++)
	{
		bird[i].x = myran->doub() * Lx;
		bird[i].y = myran->doub() * Ly;
		double theta = myran->doub() * 2 * PI;
		bird[i].vx = bird[i].vx0 = cos(theta);
		bird[i].vy = bird[i].vy0 = sin(theta);
		bird[i].next = nullptr;
	}
	return bird;
}

Node * Node::ini_move_left(Ran *myran)
{
	N = int(rho_0 * Lx * Ly);
	Node *bird = new Node[N];
	for (int i = 0; i < N; i++)
	{
		bird[i].x = myran->doub() * Lx;
		bird[i].y = myran->doub() * Ly;
		bird[i].vx = bird[i].vx0 = 1;
		bird[i].vy = bird[i].vy0 = 0;
		bird[i].next = nullptr;
	}
	return bird;
}

Node * Node::ini_copy_snap(double _eta, double _eps, double _rho_0, double _Lx, double _Ly, unsigned long long _seed, int _t)
{
	int nrows;
	int ncols;
	// check the compatibility of size of input snapshot
	if (int(Lx) % int(_Lx) == 0 && int(Ly) % int(_Ly) == 0 && rho_0 == _rho_0)
	{
		ncols = int(Lx) / int(_Lx);
		nrows = int(Ly) / int(_Ly);
	}
	else
	{
		cout << "Error, the size of input snapshot is not right" << endl;
		exit(0);
	}
	
	// read snapshot
	int _N = int(_rho_0 * _Lx * _Ly);
	float *x = new float[_N];
	float *y = new float[_N];
	float *vx = new float[_N];
	float *vy = new float[_N];

	char seed_str[30];
	num_to_str(_seed, seed_str);
	char para[100];
	snprintf(para, 100, "%g_%g_%g_%d_%d_%s_%08d", _eta, _eps, _rho_0, int(_Lx), int(_Ly), seed_str, _t);
	char infile[100];
#ifdef _MSC_VER
	snprintf(infile, 100, "snap\\s_%s.bin", para);
#else
	snprintf(infile, 100, "snap/s_%s.bin", para);
#endif
	ifstream fin(infile, ios::binary);
	if (!fin.is_open())
	{
		cout << "Error, failed to open " << infile << endl;
		exit(0);
	}
	else
	{
		float *buff = new float[3 * _N];
		fin.read((char*)&buff[0], sizeof(float) * _N * 3);
		fin.close();
		for (int i = 0; i < _N; i++)
		{
			x[i] = buff[i * 3];
			y[i] = buff[i * 3 + 1];
			float theta = buff[i * 3 + 2];
			vx[i] = cos(theta);
			vy[i] = sin(theta);
		}
		delete[] buff;
		buff = nullptr;
	}

	// initilaizing
	N = int(rho_0 * Lx * Ly);
	Node *bird = new Node[N];

	for (int row = 0; row < nrows; row++)
	{
		double dy = row * _Ly;
		for (int col = 0; col < ncols; col++)
		{
			double dx = col * _Lx;
			for (int i = 0; i < _N; i++)
			{
				int j = i + (col + row * ncols) * _N;
				bird[j].x = x[i] + dx;
				bird[j].y = y[i] + dy;
				bird[j].vx = bird[j].vx0 = vx[i];
				bird[j].vy = bird[j].vy0 = vy[i];
				bird[j].next = nullptr;
			}
		}
	}
	delete[] x;
	delete[] y;
	delete[] vx;
	delete[] vy;
	x = y = vx = vy = nullptr;
	cout << "load snapshot successfully" << endl;
	return bird;
}
