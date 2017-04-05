#include "run.h"

using namespace std;

int main(int argc, char* argv[])
{
	// set cmdline arguments
	cmdline::parser cmd;
	cmd.add<double>("eta", '\0', "noise strength", true);
	cmd.add<double>("eps", '\0', "disorder strength", false, 0);
	cmd.add<double>("Lx", 'L', "system length in x direction", true);
	cmd.add<double>("Ly", '\0', "system length in y direction", false);
	cmd.add<int>("nstep", 'n', "total steps to run", true);
	cmd.add<unsigned long long int>("seed", 's', "seed of random number", false, 1);
	cmd.add<double>("rho0", '\0', "density", false);
	cmd.add<string>("file", 'f', "filename of the input snapshot", false, "");
	cmd.add<int>("idx_frame", '\0', "index of frame to read", false, -1);
	cmd.add<string>("snap_mode", '\0', "mode to output snapshots", false, "one", cmdline::oneof<string>("one", "mult", "none"));
	cmd.add<int>("snap_dt", '\0', "time interval between two snaps", false, 2000);
	cmd.add<string>("ini_mode", 'i', "mode to initialization", false, "rand", cmdline::oneof<string>("rand", "left"));
	cmd.parse_check(argc, argv);

	//get parameters from cmdline
	double eta = cmd.get<double>("eta");
	double epsilon = cmd.get<double>("eps");
	int nStep = cmd.get<int>("nstep");
	unsigned seed = cmd.get<unsigned long long int>("seed");
	Node::Lx = cmd.get<double>("Lx");
	Node::Ly = cmd.exist("Ly") ? cmd.get<double>("Ly") : Node::Lx;

	//initialize random number generator
	Ran *myran = new Ran(seed);

	//initialize cell list; calculate total number of cells.
	Grid *cell = Grid::ini(Node::Lx, Node::Ly);

	//initial location of birds
	Node * bird = nullptr;
	ini_birds(&bird, myran, cmd);

	//set output
	Output out(eta, epsilon, Node::rho_0, Node::Lx, Node::Ly,
		Node::N, seed, nStep, 100000, 100,
		cmd.get<int>("snap_dt"), cmd.get<string>("snap_mode"));

	//link birds to cell list
	Grid::link_nodes(cell, bird);

	cout << "eta = " << eta << endl;
	cout << "epsilon = " << epsilon << endl;
	cout << "rho_0 = " << Node::rho_0 << endl;
	cout << "Lx = " << Node::Lx << endl;
	cout << "Ly = " << Node::Ly << endl;
	cout << "seed = " << seed << endl;
	cout << "tot steps = " << nStep << endl;
	cout << "tot cells = " << Grid::mm << endl;

	//run
	double *disorder = nullptr;
	ini_rand_torques(&disorder, Grid::mm, epsilon, seed);
	run(bird, cell, myran, nStep, eta, epsilon, disorder, out);

}
