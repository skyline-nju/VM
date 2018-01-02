#include "io_data.h"
#include <iomanip>

_Writer::_Writer(const cmdline::parser & cmd): idx_frame(0) {
  Lx = cmd.get<double>("Lx");
  Ly = cmd.exist("Ly") ? cmd.get<double>("Ly") : Lx;
  rho0 = cmd.get<double>("rho0");
  nBird = int(Lx * Ly * rho0);
  eta = cmd.get<double>("eta");
  eps = cmd.get<double>("eps");
  nstep = cmd.get<int>("nstep");
  if (cmd.exist("ext_torque")) {
    flag_ext_torque = true;
    ext_torque = cmd.get<double>("ext_torque");
  } else {
    ext_torque = 0;
    flag_ext_torque = false;
  }
  seed = cmd.get<unsigned long long>("seed");
}

OrderParaWriter::OrderParaWriter(const cmdline::parser & cmd) :
                                 _Writer(cmd) {
  mkdir("phi");
  char filename[100];
  if (!flag_ext_torque) {
    snprintf(filename, 100, "phi%s%g_%g_%g_%llu.dat",
      delimiter.c_str(), Lx, eta, eps, seed);
  } else {
    snprintf(filename, 100, "phi%s%g_%g_%g_%llu_%g.dat",
      delimiter.c_str(), Lx, eta, eps, seed, ext_torque);
  }
  fout.open(filename);
  set_frames(cmd);
  std::cout << "Order parameter: " << filename << "\n";
}

void OrderParaWriter::set_frames(const cmdline::parser &cmd) {
  int dn = 100;
  int i = dn;
  while (i <= nstep) {
    frames.push_back(i);
    i += dn;
  }
}

void OrderParaWriter::write(int i, const VM *birds) {
  if (!frames.empty() && i == frames[idx_frame]) {
    idx_frame++;
    double vx_m, vy_m;
    birds->get_v_mean(vx_m, vy_m);
    double phi = std::sqrt(vx_m * vx_m + vy_m * vy_m);
    double theta = std::atan2(vy_m, vx_m);
    fout << i << "\t" << std::setprecision(8) << phi 
         << "\t" << theta << "\n";
#ifdef _MSC_VER
    std::cout << i << "\t" << std::setprecision(8) << phi << "\n";
#endif
  }
}

LogWriter::LogWriter(const cmdline::parser &cmd): _Writer(cmd) {
  mkdir("log");
  char filename[100];
  if (!flag_ext_torque) {
    snprintf(filename, 100, "log%s%g_%g_%g_%llu.dat",
      delimiter.c_str(), Lx, eta, eps, seed);
  } else {
    snprintf(filename, 100, "log%s%g_%g_%g_%llu_%g.dat",
      delimiter.c_str(), Lx, eta, eps, seed, ext_torque);
  }
  fout.open(filename);
  set_frames(cmd);

  t_start = std::chrono::system_clock::now();
  std::time_t start_time = std::chrono::system_clock::to_time_t(t_start);
  fout << "Started simulation at " << std::ctime(&start_time) << "\n";
  fout << "-------- Parameters --------\n";
  fout << "Particle number: " << Lx * Ly * rho0 << "\n";
  fout << "Density: " << rho0 << "\n";
  fout << "eta: " << eta << "\n";
  fout << "eps: " << eps << "\n";
  fout << "Lx: " << Lx << "\n";
  fout << "Ly: " << Ly << "\n";
  fout << "seed: " << seed << "\n";
  fout << "Total time steps: " << nstep << "\n";
  fout << "\n";
  fout << "-------- Mode --------\n";
  if (cmd.exist("metric_free"))
    fout << "Metric-free: ON\n";
  else
    fout << "Metric-free: OFF\n";
  if (cmd.exist("dt"))
    fout << "Noise mode: continuos, dt = " << cmd.get<double>("dt") << "\n";
  else if (cmd.exist("vec_noise"))
    fout << "Noise mode: vectorial\n";
  else
    fout << "Noise mode: scalar\n";
  if (cmd.exist("defect_sep"))
    fout << "Initial configuration: defect pair with separation "
    << cmd.get<double>("defect_sep") << " in mode "
    << cmd.get<int>("defect_mode") << "\n";
  else
    fout << "Initial configuration: random\n";
  if (cmd.exist("cg_on"))
    fout << "Coarse-grained snapshots: ON\n";
  else
    fout << "Coarse-grained snapshots: OFF\n";
  fout << "\n";
  fout << "-------- RUN --------\n";
  fout << "time step\telapsed time\n";
}

void LogWriter::set_frames(const cmdline::parser & cmd) {
  int dn = cmd.get<int>("log_dt");
  int i = dn;
  while (i <= nstep) {
    frames.push_back(i);
    i += dn;
  }
}

void LogWriter::write(int i, const VM * bird) {
  if (!frames.empty() && i == frames[idx_frame]) {
    idx_frame++;
    auto t_now = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = t_now - t_start;
    double dt = elapsed_seconds.count();
    int hour = int(dt / 3600);
    int min = int((dt - hour * 3600) / 60);
    int sec = dt - hour * 3600 - min * 60;
    fout << i << "\t" << hour << ":" << min << ":" << sec << "\n";
    if (i == nstep) {
      std::time_t end_time = std::chrono::system_clock::to_time_t(t_now);
      fout << "Finished simulation at " << std::ctime(&end_time) << "\n";
    }
  }
}

CoarseGrainSnapWriter::CoarseGrainSnapWriter(const cmdline::parser &cmd):
                                             _Writer(cmd) {
  mkdir("snap");
  char filename[100];
  l = cmd.get<double>("cg_l");
  ncols = int(Lx / l);
  int nrows = int(Ly / l);
  ncells = ncols * nrows;
  std::string prefix = "snap" + delimiter;
  dt = cmd.get<double>("dt");
  if (cmd.exist("dt")) {
    prefix += "CHff";
    dt = cmd.get<double>("dt");
    is_continous = true;
  } else {
    prefix += "cHff";
    dt = 1;
    is_continous = false;
  }

  if (!flag_ext_torque) {
    snprintf(filename, 100, "%s_%g_%g_%g_%g_%d_%d_%d_%llu.bin",
      prefix.c_str(), eta, eps, Lx, Ly, ncols, nrows, nBird, seed);
  } else {
    snprintf(filename, 100, "%s_%g_%g_%g_%g_%d_%d_%d_%llu_%g.bin",
      prefix.c_str(), eta, eps, Lx, Ly, ncols, nrows, nBird, seed, ext_torque);
  }
  fout.open(filename, std::ios::binary);
  set_frames(cmd);
  std::cout << "snapshot: ON\t " << filename << "\n";
}

void CoarseGrainSnapWriter::set_frames(const cmdline::parser &cmd) {
  if (cmd.exist("cg_dt")) {
    int dn = cmd.get<int>("cg_dt");
    int i = dn;
    while (i <= nstep) {
      frames.push_back(i);
      i += dn;
    }
  } else {
    int dn;
    int i = 1;
    while (i <= nstep) {
      frames.push_back(i);
      if (i < 20)
        dn = 1;
      else if (i < 80)
        dn = 2;
      else if (i < 160)
        dn = 4;
      else if (i < 400)
        dn = 8;
      else if (i < 1000)
        dn = 16;
      else
        dn = 32;
      i += dn;
    }
    std::cout << "i_last = " << i << "\n";
  }
}

void CoarseGrainSnapWriter::write(int i, const VM * birds) {
  if ((!frames.empty() && i == frames[idx_frame]) || i == 0) {
    if (i > 0)
      idx_frame++;
    unsigned short *n_cg = new unsigned short[ncells];
    float *vx_cg = new float[ncells];
    float *vy_cg = new float[ncells];
    double vx_m, vy_m;
    coarse_grain(n_cg, vx_cg, vy_cg, vx_m, vy_m, l, ncols, ncells, birds);
    if (is_continous) {
      double t = dt * i;
      fout.write((char*)&t, sizeof(double));
    } else {
      fout.write((char*)&i, sizeof(int));
    }
    fout.write((char*)&vx_m, sizeof(double));
    fout.write((char*)&vy_m, sizeof(double));
    fout.write((char *)n_cg, sizeof(unsigned short) * ncells);
    fout.write((char *)vx_cg, sizeof(float) * ncells);
    fout.write((char *)vy_cg, sizeof(float) * ncells);

    delete[] n_cg;
    delete[] vx_cg;
    delete[] vy_cg;
#ifdef _MSC_VER
    std::cout << "output coarse-grained snapshot at step " << i << "\n";
#endif
  }
}


