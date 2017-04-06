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
	cmd.add<unsigned long long int>(
		"seed", 's', "seed of random number", false, 1);
	cmd.add<double>("rho0", '\0', "density", false);
	cmd.add<string>("file", 'f', "input file", false, "");
	cmd.add<int>("idx_frame", '\0', "which frame to read", false, -1);
	cmd.add<string>("snap_mode", '\0', "mode to output snapshots", false, 
					"one", cmdline::oneof<string>("one", "mult", "none"));
	cmd.add<int>(
		"snap_dt", '\0', "time interval to output snap", false, 2000);
	cmd.add<int>(
		"log_dt", '\0', "time interval to log", false, 100000);
	cmd.add<int>(
		"phi_dt", '\0', "time interval to calculate phi", false, 100);
	cmd.add<string>("ini_mode", 'i', "initializing mode", false, 
		"rand", cmdline::oneof<string>("rand", "left"));
	cmd.add<int>("cg_dt", '\0', "time interval to coarse grain", false, 0);
	cmd.add<int>("cg_ncol", '\0', "number of cols for coarse grain", false);
	cmd.add<int>("cg_nrow", '\0', "number of rows for coarse grain", false);
	cmd.add<double>("cg_lx", '\0', "Box size in x for coarse grain", false);
	cmd.add<double>("cg_ly", '\0', "Box size in y for coarse grain", false);
	cmd.add<string>(
		"cg_format", '\0', "file format for coarse grain", false, "iff");
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
	Node *bird = nullptr;
	ini_birds(&bird, myran, cmd);

	//set output
	Output out(Node::rho_0, Node::Ly, Node::N, cmd);

	//link birds to cell list
	Grid::link_nodes(cell, bird);

	//initialize disorder
	double *disorder = nullptr;
	ini_rand_torques(&disorder, Grid::mm, epsilon, seed);
	//run
	run(bird, cell, myran, nStep, eta, epsilon, disorder, out);

}
