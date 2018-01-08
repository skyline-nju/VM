#include "cmdline.h"
#include "io_data.h"

int main(int argc, char* argv[]) {
  //set parameters
  cmdline::parser cmd;
  cmd.add<double>("eta", '\0', "noise strength", true);
  cmd.add<double>("eps", '\0', "disorder strength", false, 0);
  cmd.add<double>("Lx", 'L', "system length in x direction", true);
  cmd.add<double>("Ly", '\0', "system length in y direction", false);
  cmd.add<double>("rho0", '\0', "particle density", false, 1);
  cmd.add<double>("v0", '\0', "fixed velocity", false, 0.5);
  cmd.add<int>("nstep", 'n', "total steps to run", true);
  cmd.add<double>("dt", '\0', "time interval between two steps", false, 1);
  cmd.add<double>("ext_torque", '\0', "extra torque", false, 0);
  cmd.add<unsigned long long int>(
    "seed", 's', "seed of random number", false, 1);
  cmd.add<double>("defect_sep", '\0', "separation between two defects", false, 0);
  cmd.add<int>("defect_mode", '\0', "defect mode", false, 0);
  cmd.add<int>("log_dt", '\0', "step interval to record log", false, 1000);
  // coarse-grained snapshots
  cmd.add("cg_on", '\0', "output coarse-grained snapshots");
  cmd.add<double>("cg_l", '\0', "boxes size for coarse graining", false, 1);
  cmd.add<int>("cg_dt", '\0', "step interval to output snapshot", false, 10000);
  // whether to use vectorial noise
  cmd.add("vec_noise", '\0', "use vectorial noise");
  // is metric-free
  cmd.add("metric_free", '\0', "is metric-free?");
  cmd.parse_check(argc, argv);

  Ran myran(cmd.get<unsigned long long>("seed"));
  VM *birds = NULL;
  if (cmd.exist("metric_free")) {
    if (cmd.exist("dt")) {
      birds = new VM_metric_free<V_conti>(cmd, myran);
      std::cout << "Continous step, dt = " << cmd.get<double>("dt") << "\n";
    } else if (cmd.exist("vec_noise")) {
      birds = new VM_metric_free<V_vectorial>(cmd, myran);
      std::cout << "Vectorial Noise, dt = " << cmd.get<double>("dt") << "\n";
    } else {
      birds = new VM_metric_free<V_scalar>(cmd, myran);
      std::cout << "Scalar Noise, dt = " << cmd.get<double>("dt") << "\n";
    }
  } else {
    if (cmd.exist("dt")) {
      birds = new VM_metric<V_conti>(cmd, myran);
      std::cout << "Continous step, dt = " << cmd.get<double>("dt") << "\n";
    } else if (cmd.exist("vec_noise")) {
      birds = new VM_metric<V_vectorial>(cmd, myran);
      std::cout << "Vectorial Noise, dt = " << cmd.get<double>("dt") << "\n";
    } else {
      birds = new VM_metric<V_scalar>(cmd, myran);
      std::cout << "Scalar Noise, dt = " << cmd.get<double>("dt") << "\n";
    }
  }

  ini_output(cmd, birds);
  int n = cmd.get<int>("nstep");
  double dt = cmd.get<double>("dt");
  for (int i = 1; i <= n; i++) {
    birds->align();
    birds->stream(dt, myran);
    output(i, birds);
  }
  delete birds;
}
