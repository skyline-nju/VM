#include "run.h"
#include "iodata.h"
#include "cmdline.h"
#include <iostream>

using namespace std;

int main(int argc, char* argv[])
{
	// define parameters
	double eta;
	double epsilon;
	int nStep;
	unsigned long long seed;
	Ran *myran = nullptr;
	double *disorder = nullptr;
	Node *bird = nullptr;
	Grid *cell = nullptr;
	int dt_snap = 2000;
	string infile;
	bool flag_ini_rand = true;
	bool flag_ini_move_left = false;
	bool flag_ini_copy_snap = false;
	bool flag_disorder = true;
	bool flag_one_snap_file = true;

#ifdef _MSC_VER
	eta = 0.2;
	epsilon = 0;
	Node::rho_0 = 1.0;
	Node::Lx = 100;
	Node::Ly = 100;
	seed = 123;
	nStep = 200000;
#else
	cmdline::parser cmd;
	cmd.add<double>("eta", '\0', "noise strength", true);
	cmd.add<double>("eps", '\0', "disorder strength", true);
	cmd.add<double>("Lx", '\0', "system length in x direction", true);
	cmd.add<double>("Ly", '\0', "system length in y direction", false);
	cmd.add<int>("nstep", '\0', "total steps to run", true);
	cmd.add<unsigned long long int>("seed", '\0', "seed of random number", false, 1);
	cmd.add<double>("rho0", '\0', "density", false);
	cmd.add<string>("file", 'f', "filename of the input snapshot", false, "");
	cmd.add<int>("dt_snap", '\0', "time interval between two snap", false, 2000);
	cmd.add("moveleft", '\0', "all move toward left at beginning");
	cmd.add("multsnap", '\0', "multiple snap files");

	cmd.parse_check(argc, argv);

	eta = cmd.get<double>("eta");
	epsilon = cmd.get<double>("eps");
	nStep = cmd.get<int>("nstep");
	seed = cmd.get<unsigned long long int>("seed");
	Node::Lx = cmd.get<double>("Lx");
	dt_snap = cmd.get<int>("dt_snap");
	if (cmd.exist("Ly"))
		Node::Ly = cmd.get<double>("Ly");
	else
		Node::Ly = Node::Lx;

	if (!cmd.exist("rho0") && !cmd.exist("file"))
	{
		cerr<<cmd.usage();
		return 0;
	}
	else if (cmd.exist("file"))
	{
		infile = cmd.get<string>("file");
		flag_ini_copy_snap = true;
	}
	else
	{
		Node::rho_0=cmd.get<double>("rho0");
		if (cmd.exist("moveleft"))
			flag_ini_move_left = true;
		else
			flag_ini_rand = true;
	}

	if (cmd.exist("multsnap"))
		flag_one_snap_file = false;
#endif

	//initialize random number generator
	myran = new Ran(seed);

	//initialize cell list; calculate total number of cells.
	cell = Grid::ini(Node::Lx, Node::Ly);

	//initialize the quenched disorder
	if (flag_disorder)
		ini_rand_torques(&disorder, Grid::mm, epsilon, myran);

	//initialize the coordination of birds
	if (flag_ini_rand)
		bird = Node::ini_rand(myran);
	else if (flag_ini_move_left)
		bird = Node::ini_move_left(myran);

	//set output
	ini_output(eta, epsilon, seed, nStep, Grid::mm, flag_one_snap_file, dt_snap);

	//link birds to cell list
	Grid::link_nodes(cell, bird);

	//run
	if (flag_disorder)
		run(bird, cell, myran, nStep, eta, epsilon, disorder);
	else
		run(bird, cell, myran, nStep, eta);
}
