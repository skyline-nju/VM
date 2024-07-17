#include "io_data.h"
#include <iomanip>

std::vector<_Writer *> writers;
// std::vector<std::ofstream> fouts;

// parameters
double Lx;
double Ly;
double rho0;
int nBird;
int nstep;
double eta;
double eps;
double ext_torque;
unsigned long long seed;
bool flag_ext_torque;
Snap_GSD_2 *gsd_snap = nullptr;

void ini_output(const cmdline::parser &cmd, const VM *birds, const std::string &prefix) {
  // char prefix[100] = ".";
  // set parameters
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

  // initialize writers
  if (cmd.exist("cg_on")) {
    writers.push_back(new CoarseGrainSnapWriter(cmd, prefix.c_str()));
    writers[0]->write(0, birds);
  }
  writers.push_back(new OrderParaWriter(cmd, prefix));
  writers.push_back(new LogWriter(cmd, prefix));

  if (cmd.exist("snap_dt")) {
    char filename[255];
    if (cmd.exist("alpha")) {
      snprintf(filename, 255, "%s%sL%g_%g_a%g_e%.3f_r%g_s%llu.gsd", prefix.c_str(), delimiter.c_str(), Lx, Ly, cmd.get<double>("alpha"), eta, rho0, seed);
    } else if (cmd.exist("dis_frac")) {
      snprintf(filename, 255, "%s%sL%g_%g_d%.4f_e%.3f_r%g_s%llu.gsd", prefix.c_str(), delimiter.c_str(), Lx, Ly, cmd.get<double>("dis_frac"), eta, rho0, seed);
    } else {
      snprintf(filename, 255, "%s%sL%g_%g_e%.3f_r%g_s%llu.gsd", prefix.c_str(), delimiter.c_str(), Lx, Ly, eta, rho0, seed);
    }
    gsd_snap = new Snap_GSD_2(filename, 0, nstep, cmd.get<int>("snap_dt"), Lx, Ly, "rand");
  }
}

void output(int i, const VM* birds) {
  for (int j = 0; j < writers.size(); j++) {
    writers[j]->write(i, birds);
  }
  if (gsd_snap) {
    gsd_snap->dump(i, birds);
  }
}

OrderParaWriter::OrderParaWriter(const cmdline::parser & cmd, const std::string& prefix) : _Writer(cmd, prefix) {
  char folder[255];
  snprintf(folder, 255, "%s%sphi%s", prefix.c_str(), delimiter.c_str(), delimiter.c_str());
  mkdir(folder);
  char filename[255];
  if (cmd.exist("alpha")) {
    snprintf(filename, 255, "%s%g_%g_%g_%g_%g_%llu.dat", folder, Lx, Ly, cmd.get<double>("alpha"), eta, rho0, seed);
  } else if (cmd.exist("dis_frac")) {
    snprintf(filename, 255, "%s%g_%g_%.4f_%.3f_%.3f_%llu.dat", folder, Lx, Ly, cmd.get<double>("dis_frac"), eta, rho0, seed);
  } else {
    snprintf(filename, 255, "%s%g_%g_%g_%g_%llu.dat", folder, Lx, Ly, eta, rho0, seed);
  }
  fout_.open(filename);
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
    fout_ << i << "\t" << std::setprecision(8) << phi 
         << "\t" << theta << std::endl;
#ifdef _MSC_VER
    std::cout << i << "\t" << std::setprecision(8) << phi << "\n";
#endif
  }
}

LogWriter::LogWriter(const cmdline::parser &cmd, const std::string& prefix): _Writer(cmd, prefix) {
  char folder[255];
  snprintf(folder, 255, "%s%slog%s", prefix.c_str(), delimiter.c_str(), delimiter.c_str());
  mkdir(folder);
  char filename[255];
  if (cmd.exist("alpha")) {
    snprintf(filename, 255, "%s%g_%g_%g_%g_%g_%llu.dat", folder, Lx, Ly, cmd.get<double>("alpha"), eta, rho0, seed);
  } else if (cmd.exist("dis_frac")) {
    snprintf(filename, 255, "%s%g_%g_%g_%g_%g_%llu.dat", folder, Lx, Ly, cmd.get<double>("dis_frac"), eta, rho0, seed);
  } else {
    snprintf(filename, 255, "%s%g_%g_%g_%g_%llu.dat", folder, Lx, Ly, eta, rho0, seed);
  }
  fout_.open(filename);
  set_frames(cmd);

  t_start = std::chrono::system_clock::now();
  std::time_t start_time = std::chrono::system_clock::to_time_t(t_start);
  fout_ << "Started simulation at " << std::ctime(&start_time) << "\n";
  fout_ << "-------- Parameters --------\n";
  fout_ << "Particle number: " << Lx * Ly * rho0 << "\n";
  fout_ << "Density: " << rho0 << "\n";
  fout_ << "eta: " << eta << "\n";
  fout_ << "eps: " << eps << "\n";
  fout_ << "Lx: " << Lx << "\n";
  fout_ << "Ly: " << Ly << "\n";
  fout_ << "seed: " << seed << "\n";
  fout_ << "Total time steps: " << nstep << "\n";
  fout_ << "\n";
  fout_ << "-------- Mode --------\n";
  if (cmd.exist("metric_free"))
    fout_ << "Metric-free: ON\n";
  else
    fout_ << "Metric-free: OFF\n";
  if (cmd.exist("dt"))
    fout_ << "Noise mode: continuos, dt = " << cmd.get<double>("dt") << "\n";
  else if (cmd.exist("vec_noise"))
    fout_ << "Noise mode: vectorial\n";
  else
    fout_ << "Noise mode: scalar\n";
  if (cmd.exist("alpha")) {
    fout_ << "Fore-arf asymmetry with alpha=" << cmd.get<double>("alpha") << "\n";
  }
  if (cmd.exist("defect_sep"))
    fout_ << "Initial configuration: defect pair with separation "
    << cmd.get<double>("defect_sep") << " in mode "
    << cmd.get<int>("defect_mode") << "\n";
  else
    fout_ << "Initial configuration: random\n";
  if (cmd.exist("cg_on"))
    fout_ << "Coarse-grained snapshots: ON\n";
  else
    fout_ << "Coarse-grained snapshots: OFF\n";
  fout_ << "\n";
  fout_ << "-------- RUN --------\n";
  fout_ << "time step\telapsed time" << std::endl;
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
    fout_ << i << "\t" << hour << ":" << min << ":" << sec << std::endl;
    if (i == nstep) {
      std::time_t end_time = std::chrono::system_clock::to_time_t(t_now);
      fout_ << "Finished simulation at " << std::ctime(&end_time) << "\n";
    }
  }
}

CoarseGrainSnapWriter::CoarseGrainSnapWriter(const cmdline::parser &cmd, const std::string& prefix):
                                             _Writer(cmd, prefix) {
  char folder[255];
  snprintf(folder, 255, "%s%ssnap%s", prefix.c_str(), delimiter.c_str(), delimiter.c_str());
  mkdir(folder);
  char filename[255];
  l = cmd.get<double>("cg_l");
  ncols = int(Lx / l);
  int nrows = int(Ly / l);
  ncells = ncols * nrows;

  snprintf(filename, 256, "%s%g_%g_%g_%d_%d_%d_%llu.bin",
    folder, eta, Lx, Ly, ncols, nrows, nBird, seed);
  
  fout_.open(filename, std::ios::binary);
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
      fout_.write((char*)&t, sizeof(double));
    } else {
      fout_.write((char*)&i, sizeof(int));
    }
    fout_.write((char*)&vx_m, sizeof(double));
    fout_.write((char*)&vy_m, sizeof(double));
    fout_.write((char *)n_cg, sizeof(unsigned short) * ncells);
    fout_.write((char *)vx_cg, sizeof(float) * ncells);
    fout_.write((char *)vy_cg, sizeof(float) * ncells);

    delete[] n_cg;
    delete[] vx_cg;
    delete[] vy_cg;
#ifdef _MSC_VER
    std::cout << "output coarse-grained snapshot at step " << i << "\n";
#endif
  }
}

void ExporterBase::set_lin_frame(int start, int n_step, int sep) {
  n_step_ = n_step;
  for (auto i = start + sep; i <= n_step_; i += sep) {
    frames_arr_.push_back(i);
  }
  frame_iter_ = frames_arr_.begin();
}

bool ExporterBase::need_export(int i_step) {
  bool flag = false;
  if (!frames_arr_.empty() && i_step == (*frame_iter_)) {
    frame_iter_++;
    flag = true;
  }
  return flag;
}

Snap_GSD_2::Snap_GSD_2(const std::string& filename,
                                 int start, int n_step, int sep,
                                 double Lx, double Ly,
                                 const std::string& open_flag)
  : ExporterBase(start, n_step, sep) {
  unsigned int version = gsd_make_version(1, 4);
  handle_ = new gsd_handle;
  if (open_flag == "rand" || open_flag == "ordered") {
    int flag = gsd_create(filename.c_str(), "cpp", "hoomd", version);
    if (flag != 0) {
      std::cout << "Error when create " << filename << "; state=" << flag << std::endl;
      exit(1);
    }
    flag = gsd_open(handle_, filename.c_str(), GSD_OPEN_READWRITE);
    if (flag != 0) {
      std::cout << "Error when open " << filename << "; state=" << flag << std::endl;
      exit(1);
    }

    float box[6] = {float(Lx), float(Ly), 1, 0, 0, 0 };
    gsd_write_chunk(handle_, "configuration/box", GSD_TYPE_FLOAT, 6, 1, 0, box);
    
    char types[] = {'A', 'B'};
    gsd_write_chunk(handle_, "particles/types", GSD_TYPE_INT8, 2, 1, 0, types);
  } else if (open_flag == "resume") {
    int flag = gsd_open(handle_, filename.c_str(), GSD_OPEN_READWRITE);
    if (flag != 0) {
      std::cout << "Error when open " << filename << "; state=" << flag << std::endl;
      exit(1);
    } else {
      std::cout << "open " << filename << std::endl;
    } 
  } else {
    std::cout << "Wrong open flag, which must be one of 'rand', 'ordered' and 'resume'!" << std::endl;
    exit(1);
  }
  half_Lx_ = Lx / 2;
  half_Ly_ = Ly / 2;
}

Snap_GSD_2::~Snap_GSD_2() {
  gsd_close(handle_);
  delete handle_;
}


uint64_t Snap_GSD_2::get_time_step() {
  uint64_t step;
  size_t n_frame = gsd_get_nframes(handle_);
  if (n_frame == 0) {
    step = sep_;
  } else {
    const gsd_index_entry* chunk = gsd_find_chunk(handle_, n_frame - 1, "configuration/step");
    if (chunk) {
      gsd_read_chunk(handle_, &step, chunk);
      step += sep_;
    } else {
      step = sep_;
    }
  }
  return step;
}


 void Snap_GSD_2::get_data_from_par(const VM* birds, float* pos,
                                    uint32_t* type_id) const {
  int N = birds->get_num_birds();
  for (int i = 0; i < N; i++) {
    double x, y, theta;
    birds->get_x(i, x, y);
    birds->get_theta(i, theta);
    size_t j3 = i * 3;
    pos[j3    ] = x - half_Lx_;
    pos[j3 + 1] = y - half_Ly_;
    pos[j3 + 2] = theta;
    if (type_id != nullptr) {
        birds->get_type(i, type_id[i]);
    }
  }
 }

  void Snap_GSD_2::dump(int i_step, const VM*birds) {
    if (need_export(i_step)) {
    uint32_t n_par = birds->get_num_birds();
    float* pos = new float[n_par * 3];
    uint32_t* type_id = new uint32_t[n_par];

    get_data_from_par(birds, pos, type_id);
    uint64_t step = get_time_step();
  
    // std::cout << "dump frame " << step << std::endl;
    gsd_write_chunk(handle_, "configuration/step", GSD_TYPE_UINT64, 1, 1, 0, &step);
    gsd_write_chunk(handle_, "particles/N", GSD_TYPE_UINT32, 1, 1, 0, &n_par);
    gsd_write_chunk(handle_, "particles/position", GSD_TYPE_FLOAT, n_par, 3, 0, pos);
    gsd_write_chunk(handle_, "particles/typeid", GSD_TYPE_UINT32, n_par, 1, 0, type_id);

    gsd_end_frame(handle_);
    delete[] pos;
    delete[] type_id;
  }
}
