#include "run.h"
#include "Par.h"
#include "exporter.h"

#include <iostream>
int main(int argc, char* argv[]) {
#ifdef _MSC_VER
  double Lx = 128;
  double Ly = 128;
  double rho0 = 0.5;
  double eta = 0.2;
  double v0 = 0.5;

  int n_step = 10000;
  int snap_interval = 1000;
  int seed = 100;

  std::string ini_mode = "resume";
  char folder[] = "D:/code/VM/AsymVM/data/";
#else
  double Lx = atof(argv[1]);
  double Ly = atof(argv[2]);
  double rho0 = atof(argv[3]);
  double eta = atof(argv[4]);
  double v0 = 0.5;
  int n_step = atoi(argv[5]);
  int snap_interval = atoi(argv[6]);
  int seed = atoi(argv[7]);
  std::string ini_mode = argv[8];
  char folder[] = "."
#endif
  int n_par = int(round(Lx * Ly * rho0));

  std::vector<Par> p_arr;
  Grid<Par> cells(Lx, Ly, 1.0);

  Ranq2 myran(seed);

  // set output
  char basename[255];
  char log_file[255];
  char gsd_file[255];
  snprintf(basename, 255, "L%g_%g_e%g_r%g_s%d", Lx, Ly, eta, rho0, seed);
  snprintf(log_file, 255, "%slog_%s.dat", folder, basename);
  snprintf(gsd_file, 255, "%s%s.gsd", folder, basename);

  int log_interval = 10000;
  exporter::LogExporter log(log_file, 0, n_step, log_interval, n_par);
  exporter::Snap_GSD_2 gsd(gsd_file, 0, n_step, snap_interval, Lx, Ly, ini_mode);

  // initialize particles
  ini_particles(p_arr, myran, cells, n_par, Lx, Ly, ini_mode, gsd);


  for (int i = 1; i <= n_step; i++) {
    run_one_step(p_arr, myran, cells, eta, Lx, Ly, v0);
    log.record(i);
    gsd.dump(i, p_arr);
  }

}
