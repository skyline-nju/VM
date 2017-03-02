#include "disorder.h"
#include "output.h"

using namespace std;

double eta;
double epsilon;
int Nstep;
int seed;

Ran *myran = NULL;
double *disorder = NULL;
Node *bird = NULL;
Grid *cell = NULL;

void ini(int argc, char ** argv)
{
#ifdef _MSC_VER
	eta = 0.2;
	epsilon = 0;
	Node::rho_0 = 1.0;
	Node::Lx = 200;
	Node::Ly = 200;
	seed = 123;
	Nstep = 200000;
#else
	if (argc < 8)
	{
		cout << "Error, need 7 arguments at least, but only " << argc << " was given" << endl;
		cout << "Argument list: " << endl;
		cout << "eta, epsilon, rho_0, Lx, Ly, seed, nStep" << endl;
		cout << "eta, epsilon, rho_0, Lx, Ly, seed, nStep, t" << endl;
		cout << "eta, epsilon, rho_0, Lx1, Ly1, seed1, nStep, Lx2, Ly2, seed2, t" << endl;
		exit(0);
	}
	eta = atof(argv[1]);
	epsilon = atof(argv[2]);
	Node::rho_0 = atof(argv[3]);
	Node::Lx = atof(argv[4]);
	Node::Ly = atof(argv[5]);
	seed = atoi(argv[6]);
	Nstep = atoi(argv[7]);
#endif

	myran = new Ran(seed);

	/* set mx, my, mm */
	cell = Grid::ini();

	/* initialize quenched disorder */
	disorder = ini_torques(myran, epsilon);

	/* initialize location and velocity of birds */
#ifdef _MSC_VER
	//bird = Node::ini_move_left(myran);
	bird = Node::ini_snap(eta, epsilon, Node::rho_0, 0.5 * Node::Lx, 0.5 * Node::Ly, seed, 100000);
#else
	if (argc == 8)
	{
		bird = Node::ini_move_left(myran);
	}
	else if (argc == 9)
	{
		int t0 = atoi(argv[8]);
		bird = Node::ini_snap(eta, epsilon, Node::rho_0, Node::Lx, Node::Ly, seed, t0);
	}
	else if (argc == 12)
	{
		double Lx2 = atof(argv[8]);
		double Ly2 = atof(argv[9]);
		int seed2 = atoi(argv[10]);
		int t2 = atoi(argv[11]);
		bird = Node::ini_snap(eta, epsilon, Node::rho_0, Lx2, Ly2, seed2, t2);
	}
	else if (argc == 15)
	{
		double eta2 = atof(argv[8]);
		double eps2 = atof(argv[9]);
		double rho2 = atof(argv[10]);
		double Lx2 = atof(argv[11]);
		double Ly2 = atof(argv[12]);
		int seed2 = atoi(argv[13]);
		int t2 = atoi(argv[14]);
		bird = Node::ini_snap(eta2, eps2, rho2, Lx2, Ly2, seed2, t2);
	}
	else
	{
		cout << "Error, improper number of argments" << endl;
	}
#endif

	//initialize outputting
	output_ini(eta, epsilon, seed, Nstep);

	//link birds to cells
	Grid::link_nodes(cell, bird);
}

void update()
{
	static double eta2PI = eta * 2 * PI;

	//updating velocity
	Grid::all_pairs(cell);

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
	noise = NULL;

	//updating cell list
	Grid::refresh(cell, bird);
}

void run(int tot_step)
{
	for (int i = 0; i < tot_step; i++)
	{
		output(bird, i);
		update();
	}
	output(bird, tot_step);
}

int main(int argc, char* argv[])
{
	ini(argc, argv);
	run(Nstep);
	return 0;
}
