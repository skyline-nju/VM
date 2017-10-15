#include "io_data.h"
#include <cmath>

using namespace std;

#ifdef _MSC_VER
string delimiter("\\");
#else
string delimiter("/");
#endif

void get_mean_velocity(const Node *bird, int nBird, double &vxm, double &vym) {
  vxm = 0;
  vym = 0;
  for (int i = 0; i < nBird; i++) {
    vxm += bird[i].vx;
    vym += bird[i].vy;
  }
  vxm /= nBird;
  vym /= nBird;
}

void gene_log_frames(vector<int> &frames, double exponent, int t_end) {
  frames.push_back(1);
  double t = 1;
  int cur_frame = 1;
  while (true) {
    t *= exponent;
    if (t >= t_end) break;
    int tmp_frame = round(t);
    if (tmp_frame > cur_frame) {
      cur_frame = tmp_frame;
      frames.push_back(cur_frame);
    }
  }
}

void gene_log_frames(vector<int> &frames, double exponent, int t_beg, int t_end) {
  int t = t_beg;
  while (t <= t_end) {
    frames.push_back(t);
    t *= exponent;
  }
}

void gene_log_frames_with_log_win(vector<int> &frames) {
  int t_mid = 25;
  int win = 1;
  frames.push_back(t_mid);
  for (int i = 1; i <= 7; i++) {
    t_mid *= 2;
    win *= 2;
    for (int j = -win / 2; j < win / 2; j++) {
      frames.push_back(t_mid + j);
    }
  }
}

void gene_lin_frames(vector<int> &frames, int dt, int t_end) {
  for (int t = 1; t <= t_end; t++) {
    if (t % dt == 0)
      frames.push_back(t);
  }
}

void gene_lin_frames(vector<int> &frames, int dt, int t_end, int t_beg) {
  for (int t = t_beg; t <= t_end; t++) {
    if (t % dt == 0)
      frames.push_back(t);
  }
}


OrderPara::OrderPara(double eta, double eps, double rho0,
                     double Lx, double Ly, 
                     unsigned long long seed, int phi_dt, ofstream &log) {
  mkdir("phi");
  char buff[100];
  snprintf(buff, 100, "phi%sp_%g_%g_%g_%g_%g_%llu.dat",
           delimiter.c_str(), eta, eps, rho0, Lx, Ly, seed);
  fout.open(buff);
  dt = phi_dt;
  log << "file for order parameters: " << buff << endl;
  log << "interval for recording phi: " << dt << endl;
}

void OrderPara::out(const Node * bird, int nBird, int step) {
  if (step % dt == 0) {
    double vxm = 0;
    double vym = 0;
    get_mean_velocity(bird, nBird, vxm, vym);
    double phi = sqrt(vxm * vxm + vym * vym);
    double theta = atan2(vym, vxm);
    fout << fixed << std::setw(16) << setprecision(10)
      << step << "\t" << phi << "\t" << theta << endl;
  }
}

OutSnapshot::OutSnapshot(double _eta, double _eps, double _rho0,
                         double _Lx, double _Ly, unsigned long long _seed,
                         int nBird, const cmdline::parser &cmd,
                         ofstream &log, bool &flag)
    :eta(_eta), eps(_eps), rho0(_rho0), Lx(_Lx), Ly(_Ly), seed(_seed) {
  string snap_mode = cmd.get<string>("snap_mode");
  out_dt = cmd.get<int>("snap_dt");
  if (snap_mode == "one") {
    flag = true;
    is_one_file = true;
    mkdir("snap_one");
    char filename[100];
    snprintf(filename, 100, "snap_one%sso_%g_%g_%g_%g_%d_%d_%llu.bin",
      delimiter.c_str(), eta, eps, Lx, Ly, nBird, out_dt, seed);
    fout.open(filename, ios::binary);
    log << "snapshot: " << filename << endl;
  } else if (snap_mode == "mult") {
    flag = true;
    is_one_file = false;
    log << "snapshot saved in multiple files\n";
  } else {
    flag = false;
    log << "no outputing snapshots\n" << endl;
  }
}
void OutSnapshot::write(const Node * bird, int nBird) {
  float *buff = new float[3 * nBird];
  for (int j = 0; j < nBird; j++) {
    buff[3 * j] = bird[j].x;
    buff[3 * j + 1] = bird[j].y;
    buff[3 * j + 2] = atan2(bird[j].vy, bird[j].vx);
  }
  fout.write((char *)buff, sizeof(float) * nBird * 3);
  delete[] buff;
  buff = nullptr;
}

void OutSnapshot::to_file(const Node *bird, int nBird, int step) {
  if (step % out_dt == 0) {
    if (is_one_file)
      write(bird, nBird);
    else {
      char filename[100];
      snprintf(filename, 100, "snap%ss_%g_%g_%g_%g_%g_%llu_%08d.bin",
               delimiter.c_str(), eta, eps, rho0, Lx, Ly, seed, step);
      fout.open(filename, ios::binary);
      write(bird, nBird);
      fout.close();
    }
  }
}

InSnapshot::InSnapshot(const string infile) {
  vector<string> str_list = split(infile, "_");
  string folder;
  if (str_list[0] == "s") {
    folder = "snap";
    is_read_block = false;
    str_to_num(str_list[4], Lx);
    str_to_num(str_list[5], Ly);
    buff = nullptr;
  } else {
    folder = "snap_one";
    is_read_block = true;
    str_to_num(str_list[5], nBird);
    str_to_num(str_list[3], Lx);
    str_to_num(str_list[4], Ly);
    buff = new float[nBird * 3];
  }
  string path(folder + delimiter + infile);
  fin.open(path.c_str(), ios::binary | ios::ate);
  if (!fin.is_open()) {
    cout << "Error, failed to open " << infile << endl;
    exit(1);
  } else
    cout << "Open snapshots: " << infile << endl;
}

InSnapshot::~InSnapshot() {
  fin.close();
  delete[]buff;
  buff = nullptr;
}

void InSnapshot::read() {
  streamsize count = fin.tellg();
  nBird = count / sizeof(float) / 3;
  buff = new float[nBird * 3];
  fin.seekg(0, ios::beg);
  fin.read((char *)buff, count);
}

void InSnapshot::read_block(int idx_frame) {
  streamsize count = fin.tellg();
  streamsize frame_size = sizeof(float) * nBird * 3;
  if (idx_frame < 0)
    fin.seekg(frame_size * idx_frame, ios::end);
  else
    fin.seekg(frame_size * idx_frame, ios::beg);
  fin.read((char *)buff, frame_size);
}

void InSnapshot::from_file(int idx_frame, double & _Lx, double & _Ly,
                           vector<float> &x, vector<float> &y,
                           vector<float> &theta) {
  if (is_read_block)
    read_block(idx_frame);
  else
    read();
  _Lx = Lx;
  _Ly = Ly;

  x.reserve(nBird);
  y.reserve(nBird);
  theta.reserve(nBird);
  for (int i = 0; i < nBird; i++) {
    x.push_back(buff[3 * i]);
    y.push_back(buff[3 * i + 1]);
    theta.push_back(buff[3 * i + 2]);
  }
}

CoarseGrain::CoarseGrain(int nBird, const cmdline::parser &cmd,
                         ofstream &log, bool &flag) {
  int dt = cmd.get<int>("cg_dt");
  double exponent = cmd.get<double>("cg_exp");
  if (dt > 0 || exponent > 0 || cmd.exist("cg_win")) {
    flag = true;
    set_box_size(cmd);
    set_output(cmd, nBird);
    log << "coarse grain: " << filename << endl;
  } else {
    flag = false;
    log << "coarse grain: None\n";
  }
}

void CoarseGrain::set_output(const cmdline::parser &cmd, int nBird) {
  format = cmd.get<string>("cg_format");
  mkdir("coarse");
  int dt = cmd.get<int>("cg_dt");
  double exponent = cmd.get<double>("cg_exp");
  if (exponent > 0 && dt > 0) {
    gene_log_frames(vec_frames, exponent, cmd.get<int>("t_equil"));
    gene_lin_frames(vec_frames, dt, cmd.get<int>("nstep"), cmd.get<int>("t_equil"));
  } else if (exponent > 0) {
    gene_log_frames(vec_frames, exponent, 25, cmd.get<int>("nstep"));
  } else if (dt > 0){
    gene_lin_frames(vec_frames, dt, cmd.get<int>("nstep"), cmd.get<int>("t_equil"));
  } else {
    gene_log_frames_with_log_win(vec_frames);
    for (auto i : vec_frames) {
      cout << i << endl;
    }
  }
  snprintf(filename, 100, "coarse%sc%s_%g_%g_%g_%g_%d_%d_%d_%llu.bin",
           delimiter.c_str(), format.c_str(), cmd.get<double>("eta"),
           cmd.get<double>("eps"), cmd.get<double>("Lx"), cmd.get<double>("Ly"),
           ncols, nrows, nBird, cmd.get<unsigned long long int>("seed"));
  fout.open(filename, ios::binary);
  idx_cur_frame = 0;
}

void CoarseGrain::set_box_size(const cmdline::parser &cmd) {
  double Lx = cmd.get<double>("Lx");
  double Ly = cmd.get<double>("Ly");
  if (cmd.exist("cg_ncol")) {
    ncols = cmd.get<int>("cg_ncol");
    if (cmd.exist("cg_nrow"))
      nrows = cmd.get<int>("cg_nrow");
    else
      nrows = ncols;
    lx = Lx / ncols;
    ly = Ly / nrows;
  } else if (cmd.exist("cg_lx")) {
    lx = cmd.get<double>("cg_lx");
    if (cmd.exist("cg_ly"))
      ly = cmd.get<double>("cg_ly");
    else
      ly = lx;
    ncols = int(Lx / lx);
    nrows = int(Ly / ly);
  } else {
    cout << "Error, need input cg_nol or cg_lx\n";
    exit(1);
  }
  ncells = ncols * nrows;
}

void CoarseGrain::write(const Node *bird, int nBird, int step) {
  if (step == vec_frames[idx_cur_frame]) {
    fout.write((char *)&step, sizeof(int));
    double vxm = 0;
    double vym = 0;
    get_mean_velocity(bird, nBird, vxm, vym);
    fout.write((char *)&vxm, sizeof(double));
    fout.write((char *)&vym, sizeof(double));

    if (format == "B") {
      unsigned char *num = new unsigned char[ncells];
      coarse_grain(bird, nBird, num, ncells, ncols, nrows, lx, ly, false);
      fout.write((char *)num, sizeof(unsigned char) * ncells);
      delete[] num;
    } else if (format == "Hff") {
      unsigned short *num = new unsigned short[ncells];
      float *vx = new float[ncells];
      float *vy = new float[ncells];
      coarse_grain(bird, nBird, num, vx, vy, ncells, ncols, nrows, lx, ly);
      fout.write((char *)num, sizeof(unsigned short) * ncells);
      fout.write((char *)vx, sizeof(float) * ncells);
      fout.write((char *)vy, sizeof(float) * ncells);
      delete[] num;
      delete[] vx;
      delete[] vy;
    }
    idx_cur_frame++;
    cout << "out put " << idx_cur_frame << "th frame, "
      << "time step = " << step << endl;
  }
}

Corr_r::Corr_r(const cmdline::parser &cmd, int nBird) :
               lBox(cmd.get<double>("lBox")),
               ncols(cmd.get<double>("Lx") / cmd.get<double>("lBox")) {
  corr2d = new SpatialCorr2d(ncols, ncols, cmd.get<double>("Lx"), cmd.get<double>("Ly"));
  circle_ave = new CircleAverage(ncols, ncols, lBox, r);
  c_rho_r = NULL;
  c_v_r = NULL;
  sk_rho = NULL;
  sk_v = NULL;
  double Lx = cmd.get<double>("Lx");
  int ncols = int(Lx / cmd.get<double>("lBox"));
  ncells = ncols * ncols;
  mkdir("corr_r");
  char file[100];
  double rho0 = nBird / (cmd.get<double>("Lx") * cmd.get<double>("Ly"));
  //frames.push_back(25);
  //frames.push_back(50);
  //frames.push_back(100);
  //frames.push_back(200);
  //frames.push_back(400);
  //frames.push_back(800);
  //frames.push_back(1600);
  //frames.push_back(3200);
  gene_log_frames_with_log_win(frames);

  iter = frames.begin();
  if (cmd.exist("Sk")) {
    snprintf(file, 100, "corr_r%scrsk_%g_%g_%g_%g_%g_%llu_%d_%d.bin",
             delimiter.c_str(), cmd.get<double>("eta"), cmd.get<double>("eps"),
             rho0, Lx, Lx / ncols, cmd.get<unsigned long long int>("seed"),
             r.size(), frames.size());
    flag_sk = true;
  } else {
    snprintf(file, 100, "corr_r%scr_%g_%g_%g_%g_%g_%llu_%d_%d.bin",
             delimiter.c_str(), cmd.get<double>("eta"), cmd.get<double>("eps"),
             rho0, Lx, Lx / ncols, cmd.get<unsigned long long int>("seed"),
             r.size(), frames.size());
    flag_sk = false;
  }
  fout.open(file, ios::binary);
  fout.write((char*)&r[0], sizeof(double) * r.size());
}

Corr_r::~Corr_r() {
  fout.close();
}

void Corr_r::instant(const Node * bird, int nBird) {
  int *num = new int[ncells];
  double *vx = new double[ncells];
  double *vy = new double[ncells];
  int nrows = ncols;
  coarse_grain(bird, nBird, num, vx, vy, ncells, ncols, nrows, lBox, lBox);

  double *c_rho = (double *)fftw_malloc(sizeof(double) * ncells);
  double *c_v = (double *)fftw_malloc(sizeof(double) * ncells);
  if (flag_sk) {
    double *sk_rho_2d = new double[ncells];
    double *sk_v_2d = new double[ncells];
    corr2d->eval(num, vx, vy, c_rho, c_v, sk_rho_2d, sk_v_2d, rho_m, vx_m, vy_m);
    circle_ave->eval(sk_rho_2d, sk_rho);
    circle_ave->eval(sk_v_2d, sk_v);
    delete[] sk_rho_2d;
    delete[] sk_v_2d;
  } else {
    corr2d->eval(num, vx, vy, c_rho, c_v, rho_m, vx_m, vy_m); 
  }
  circle_ave->eval(c_rho, c_rho_r);
  circle_ave->eval(c_v, c_v_r);
  delete[] num;
  delete[] vx;
  delete[] vy;
  fftw_free(c_rho);
  fftw_free(c_v);
}

void Corr_r::output(const Node * bird, int nBird, int t) {
  if (t == *iter) {
    c_rho_r = new double[r.size()];
    c_v_r = new double[r.size()];
    if (flag_sk) {
      sk_rho = new double[r.size()];
      sk_v = new double[r.size()];
    }
    instant(bird, nBird);
    fout.write((char *)&t, sizeof(int));
    fout.write((char *)&rho_m, sizeof(double));
    fout.write((char *)&vx_m, sizeof(double));
    fout.write((char *)&vy_m, sizeof(double));
    fout.write((char *)c_rho_r, sizeof(double) * r.size());
    fout.write((char *)c_v_r, sizeof(double) * r.size());
    if (flag_sk) {
      fout.write((char *)sk_rho, sizeof(double) * r.size());
      fout.write((char *)sk_v, sizeof(double) * r.size());
      delete[]sk_rho;
      delete[]sk_v;
    }
    delete[] c_rho_r;
    delete[] c_v_r;
    ++iter;
  }
}

Output::Output(double rho0, double Ly, int nBird,
               const cmdline::parser &cmd) {
  using std::placeholders::_1;
  using std::placeholders::_2;
  using std::placeholders::_3;

  double eta = cmd.get<double>("eta");
  double eps = cmd.get<double>("eps");
  double Lx = cmd.get<double>("Lx");
  unsigned long long seed = cmd.get<unsigned long long>("seed");

  time(&beg_time);
  ini_fout(eta, eps, rho0, Lx, Ly, seed);
  interval = cmd.get<int>("log_dt");
  fout << "Started at " << asctime(localtime(&beg_time));
  fout << "-------- Parameters --------\n";
  fout << "Particle number: " << nBird << endl;
  fout << "density: " << rho0 << endl;
  fout << "eta: " << eta << endl;
  fout << "epsilon: " << eps << endl;
  fout << "Lx: " << Lx << endl;
  fout << "Ly: " << Ly << endl;
  fout << "seed: " << seed << endl;
  fout << "total time steps: " << cmd.get<int>("nstep") << endl;  
  if (cmd.exist("file"))
    fout << "ini mode: " << cmd.get<string>("file") << endl;
  else
    fout << "ini mode: " << cmd.get<string>("ini_mode") << endl;
  fout << endl;
  fout << "-------- Output setting --------\n";
  phi = new OrderPara(
      eta, eps, rho0, Lx, Ly, seed, cmd.get<int>("phi_dt"), fout);
  fout_vec.push_back(bind(&OrderPara::out, phi, _1, _2, _3));

  bool flag;
  snap = new OutSnapshot(
      eta, eps, rho0, Lx, Ly, seed, nBird, cmd, fout, flag);
  if (flag)
    fout_vec.push_back(bind(&OutSnapshot::to_file, snap, _1, _2, _3));
  cg = new CoarseGrain(nBird, cmd, fout, flag);
  if (flag)
    fout_vec.push_back(bind(&CoarseGrain::write, cg, _1, _2, _3));
  if (cmd.get<double>("lBox") > 0) {
    cr = new Corr_r(cmd, nBird);
    fout_vec.push_back(bind(&Corr_r::output, cr, _1, _2, _3));
  }

  fout << endl;
  fout << "-------- Run --------\n";
  fout << "time step\telapsed time\n";
  cout << "size of int: " << sizeof(int) << endl;
  cout << "size of double: " << sizeof(double) << endl;
}

Output::~Output() {
  if (phi) delete phi;
  if (snap) delete snap;
  if (cg) delete cg;
  time_t end_time = time(nullptr);
  fout << "Finished at " << asctime(localtime(&end_time));
  fout.close();
}

void Output::ini_fout(double eta, double eps, double rho0,
                      double Lx, double Ly, unsigned long long seed) {
  mkdir("log");
  char logfile[100];
  snprintf(logfile, 100, "log%sl_%g_%g_%g_%g_%g_%llu.dat",
           delimiter.c_str(), eta, eps, rho0, Lx, Ly, seed);
  fout.open(logfile);
}

void Output::out(const Node * bird, int nBird, int step) {
  for (auto f : fout_vec) {
    f(bird, nBird, step);
  }
  if (step % interval == 0) {
    double dt = difftime(time(nullptr), beg_time);
    int hour = int(dt / 3600);
    int min = int((dt - hour * 3600) / 60);
    double sec = dt - hour * 3600 - min * 60;
    fout << step << "\t" << hour << ":" << min << ":" << sec << endl;
  }
}
