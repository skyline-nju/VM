#include "run.h"
using namespace std;

double eta;
double epsilon;
int nStep;
unsigned long long seed;
Ran *myran = nullptr;
double *disorder = nullptr;
Node *bird = nullptr;
Grid *cell = nullptr;

void ini_rand_torques(int n)
{
	double d = 1.0 / (n - 1);
	disorder = new double[n];
	for (int i = 0; i < n; i++)
		disorder[i] = (-0.5 + i *d) * epsilon * 2.0 * PI;
	shuffle(disorder, n, myran);
}

void cout_error(int n)
{
	cout << "Error, wrong arguments!" << endl;
	cout << "Valid argments: " << endl;
	cout << "\t1)eta, epsilon, rho_0, Lx, Ly, seed, nStep" << endl;
	cout << "\t2)eta, epsilon, rho_0, Lx, Ly, seed, nStep, t" << endl;
	cout << "\t3)eta, epsilon, rho_0, Lx1, Ly1, seed1, nStep, Lx2, Ly2, seed2, t" << endl;
	exit(0);
}

void ini(int argc, char ** argv)
{
	if (argc == 1)
	{
#ifdef _MSC_VER
		eta = 0.2;
		epsilon = 0;
		Node::rho_0 = 1.0;
		Node::Lx = 100;
		Node::Ly = 100;
		seed = 123;
		nStep = 200000;
#else
		cout_error(argc);
#endif
	}
	else if (argc >= 8)
	{
		str_to_num(argv[1], eta);
		str_to_num(argv[2], epsilon);
		str_to_num(argv[3], Node::rho_0);
		str_to_num(argv[4], Node::Lx);
		str_to_num(argv[5], Node::Ly);
		str_to_num(argv[6], seed);
		str_to_num(argv[7], nStep);
	}
	else
		cout_error(argc);

	myran = new Ran(seed);
	cell = Grid::ini();
	ini_rand_torques(Grid::mm);

	//initialize coordination of birds
	switch (argc)
	{
	case 1:
		bird = Node::ini_rand(myran);
		break;
	case 8:
		//bird = Node::ini_rand(myran);
		bird = Node::ini_move_left(myran);
		break;
	case 9:
		bird = Node::ini_copy_snap(eta, epsilon, Node::rho_0, Node::Lx, Node::Ly, seed, atoi(argv[8]));
		break;
	case 12:
	{
		double Lx2 = atof(argv[8]);
		double Ly2 = atof(argv[9]);
		unsigned long long seed2;
		str_to_num(argv[10], seed2);
		int t2 = atoi(argv[11]);
		bird = Node::ini_copy_snap(eta, epsilon, Node::rho_0, Lx2, Ly2, seed2, t2);
		break;
	}
	default:
		cout_error(argc);
		break;
	}

	//initialize outputting
	output_ini(eta, epsilon, seed, nStep, Grid::mm);

	//link birds to cells
	Grid::link_nodes(cell, bird);
}

void update_coor()
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

void update()
{

	//updating velocity of birds
	Grid::all_pairs(cell);

	//updating coordinate of birds
	update_coor();

	//updating cell list
	Grid::refresh(cell, bird);
}

void run()
{
	for (int i = 1; i <= nStep; i++)
	{
		update();
		output(bird, i);
	}
}

