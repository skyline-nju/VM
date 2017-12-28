#include <fstream>
#include <ctime>
#include <chrono>
#include "metric_free_vm.h"
#include "comn.h"
#include "cmdline.h"

int main(int argc, char* argv[]) {
  // set parameters
  cmdline::parser cmd;
  cmd.add<double>("eta", '\0', "noise strength", true);
  cmd.add<double>("eps", '\0', "disorder strength", false, 0);
  cmd.add<double>("Lx", 'L', "system length in x direction", true);
  cmd.add<double>("Ly", '\0', "system length in y direction", false);
  cmd.add<double>("rho0", '\0', "particle density", false, 1);
  cmd.add<int>("nstep", 'n', "total steps to run", true);
  cmd.add<unsigned long long int>(
    "seed", 's', "seed of random number", false, 1);
  cmd.add<double>("defect_sep", '\0', "separation between two defects", false, 0);
  cmd.add<int>("defect_mode", '\0', "defect mode", false, 0);
  cmd.add("output_cg", '\0', "output coarse-grained snapshots");
  cmd.add("vectorial_noise", '\0', "use vectorial noise");
  cmd.parse_check(argc, argv);

  ini(cmd);
  run(cmd.get<int>("nstep"));
  end();
}
