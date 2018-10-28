#include "run.h"

using namespace std;

int main(int argc, char* argv[]) {
  // set cmdline arguments
  cmdline::parser cmd;
  cmd.add<double>("eta", '\0', "noise strength", true);
  cmd.add<double>("eps", '\0', "disorder strength", false, 0);
  cmd.add<double>("Lx", 'L', "system length in x direction", true);
  cmd.add<double>("Ly", '\0', "system length in y direction", false);
  cmd.add<int>("nstep", 'n', "total steps to run", true);
  cmd.add<unsigned long long>("seed", 's', "seed of random number", false, 1);
  cmd.add<double>("rho0", '\0', "density", false);
  cmd.add<string>("ini_mode", '\0', "initial mode", false, "rand",
                  cmdline::oneof<string>("left", "rand"));
  cmd.add<string>("output", 'o', "path to output data", false);
  cmd.add<int>("log_dt", '\0', "time inteval to write log", false, 0);
  cmd.add<int>("order_dt", '\0', "time inteval to output order parameters", false, 100);
  cmd.add<int>("cg_dt", '\0', "time inteval to output coarse-grained snapshots",
               false, 0);
  cmd.add<int>("cg_beg", '\0', "first frame of coarse-grained snapshots",
               false, 0);
  cmd.add<int>("cg_end", '\0', "last frame of coarse-grained snapshots",
               false, 0);
  cmd.add<int>("corr_dt", '\0', "time inteval to cal correlation functions",
               false, 0);
  cmd.add<int>("corr_beg", '\0', "first frame to cal correlation functions",
               false, 0);
  cmd.add<int>("cg_ncols", '\0', "num of cols to coarse grain", false, 512);
  cmd.add<int>("cg_nrows", '\0', "num of rows to coarse grain", false, 512);

  cmd.parse_check(argc, argv);

  //get parameters from cmdline
  double eta = cmd.get<double>("eta");
  double epsilon = cmd.get<double>("eps");
  int n_step = cmd.get<int>("nstep");
  unsigned long long seed = cmd.get<unsigned long long int>("seed");
  Par::Lx = cmd.get<double>("Lx");
  Par::Ly = cmd.exist("Ly") ? cmd.get<double>("Ly") : Par::Lx;

  //initialize random number generator
  Ran *myran = new Ran(seed);

  //initialize cell list; calculate total number of cells.
  Grid *cell = Grid::ini(Par::Lx, Par::Ly);

  //initial location of birds
  Par *bird = nullptr;
  int n_par;
  ini_birds(&bird, n_par, myran, cmd);

  //link birds to cell list
  Grid::link_nodes(cell, bird);

  //initialize disorder
  double *disorder = nullptr;
  ini_rand_torques(&disorder, Grid::mm, epsilon, seed);

  // initialize otuput
  ini_output(cmd);

  //run
  run(bird, n_par, cell, myran, n_step, eta, disorder, false);

  finish_simulation();
}
