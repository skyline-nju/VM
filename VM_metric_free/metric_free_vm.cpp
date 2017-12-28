#include "metric_free_vm.h"

double Lx;
double Ly;
double eta;
double eps;
int nStep;
double rho0;
double v0 = 0.5;
bool is_scalar_noise;

Ran *myran = NULL;
PDT *DT = NULL;
std::vector<Pair_P_I> pos_idx_pair;
std::vector<Bird_wo_pos> birds;

// quenched disorder
double *random_torque = NULL;

// time
std::chrono::time_point<std::chrono::system_clock> t_start;

// order parameters
std::ofstream fout_phi;

// snapshots
std::ofstream fout_snap;
double box_size_cg = 1; // box size for coarse-graining
int n_boxes;
unsigned short *n_cg = NULL;
float * vx_cg = NULL;
float * vy_cg = NULL;
bool flag_cg;

void ini_birds_rand(unsigned int nBird, Ran *myran) {
  birds.reserve(nBird);
  pos_idx_pair.reserve(nBird);
  for (unsigned int i = 0; i < nBird; i++) {
    double x = myran->doub() * Lx;
    double y = myran->doub() * Ly;
    pos_idx_pair.push_back(std::make_pair(Point(x, y), i));
    double vx, vy;
    myran->circle_point_picking(vx, vy);
    birds.emplace_back(vx, vy);
  }
}

// create a pair of defects with separation 2a
void create_defect_pair(double a, unsigned int nBird, double rot_angle, Ran *myran) {
  double x0 = 0.5 * Lx;
  double y0 = 0.5 * Ly;
  birds.reserve(nBird);
  pos_idx_pair.reserve(nBird);
  for (unsigned int i = 0; i < nBird; i++) {
    double x = myran->doub() * Lx;
    double y = myran->doub() * Ly;
    double X = x - x0;
    double Y = y - y0;
    double phi1 = std::atan2(Y, X - a);
    double phi2 = std::atan2(Y, X + a);
    double phi = phi1 - phi2 + rot_angle;
    pos_idx_pair.push_back(std::make_pair(Point(x, y), i));
    birds.emplace_back(std::cos(phi), std::sin(phi));
  }
}

void ini(cmdline::parser &cmd) {
  // set parameters
  Lx = cmd.get<double>("Lx");
  if (cmd.exist("Ly")) {
    Ly = cmd.get<double>("Ly");
  } else {
    Ly = Lx;
  }
  eta = cmd.get<double>("eta");
  eps = cmd.get<double>("eps");
  unsigned long long seed = cmd.get<unsigned long long>("seed");
  rho0 = cmd.get<double>("rho0");
  unsigned int nBird = int(rho0 * Lx * Ly);
  t_start = std::chrono::system_clock::now();
  std::time_t start_time = std::chrono::system_clock::to_time_t(t_start);
  std::cout << "Start computation at " << std::ctime(&start_time);
  std::cout << "number of birds: " << nBird << "\n";
  std::cout << "system size: " << Lx << ", " << Ly << "\n";
  std::cout << "eta = " << eta << "\n";
  if (cmd.exist("vectorial_noise")) {
    is_scalar_noise = false;
    std::cout << "Vectorial Noise\n";
  } else {
    is_scalar_noise = true;
    std::cout << "Scalar Noise\n";
  }
  std::cout << "epsilon = " << eps << "\n";

  // initialize birds' positions and velocities
  double defect_sep = cmd.get<double>("defect_sep");
  if (defect_sep > 0) {
    int defect_mode = cmd.get<int>("defect_mode");
    double rot_angle;
    switch (defect_mode) {
    case 0:
      rot_angle = 0;
      break;
    case 1:
      rot_angle = PI * 0.5;
      break;
    case 2:
      rot_angle = PI;
      break;
    case 3:
      rot_angle = PI * 1.5;
      break;
    default:
      break;
    }
    seed = seed * 10 + defect_mode;
    if (!is_scalar_noise) {
      seed += 4;
    }
    myran = new Ran(seed);
    create_defect_pair(defect_sep, nBird, defect_mode, myran);
    std::cout << "create a pair of defects with separation " << defect_sep
      << " in mode " << defect_mode << "\n";
  } else {
    seed = seed * 10;
    if (is_scalar_noise) {
      seed += 8;
    } else {
      seed += 9;
    }
    myran = new Ran(seed);
    ini_birds_rand(nBird, myran);
    std::cout << "initialize randomly\n";
  }

  // initialize quenched disorder
  if (eps > 0) {
    int mm = int(Lx * Ly);
    random_torque = new double[mm];
    double d = 1.0 / (mm - 1);
    for (int i = 0; i < mm; i++)
      random_torque[i] = (-0.5 + i * d) * eps * 2 * PI;
    myran->shuffle(random_torque, mm);
  }

  // order parameter
  mkdir("phi");
  char filename[100];
  snprintf(filename, 100, "phi%s%g_%g_%g_%llu.dat",
    delimiter.c_str(), Lx, eta, eps, seed);
  fout_phi.open(filename);
  std::cout << "order parameter: " << filename << "\n";

  // output snapshots
  if (cmd.exist("output_cg")) {
    flag_cg = true;
    mkdir("snap");
    char fsnap[100];
    int ncols = int(Lx / box_size_cg);
    int nrows = int(Ly / box_size_cg);
    snprintf(fsnap, 100, "snap%scHff_%g_%g_%g_%g_%d_%d_%d_%llu.bin",
      delimiter.c_str(), eta, eps, Lx, Ly, ncols, nrows, nBird, seed);
    fout_snap.open(fsnap, std::ios::binary);
    n_boxes = ncols * nrows;
    n_cg = new unsigned short[n_boxes];
    vx_cg = new float[n_boxes];
    vy_cg = new float[n_boxes];
    std::cout << "snapshot: ON\t" << fsnap << "\n";
  } else {
    flag_cg = false;
    std::cout << "snapshot: OFF\n";
  }

  Bird_wo_pos::set_para(Lx, Ly, v0);

  // initialize triangulation
  PDT::Iso_rectangle domain(0, 0, Lx, Ly);
  DT = new PDT(domain);
}

void run(int n) {
  if (flag_cg) {
    output_coarse_grained_snap(0, n_cg, vx_cg, vy_cg, n_boxes, box_size_cg,
      fout_snap, pos_idx_pair, birds);
  }
  for (int i = 1; i <= n; i++) {
    align(*DT, pos_idx_pair, birds);
    move_forward(eta, *myran, is_scalar_noise, random_torque, pos_idx_pair, birds);
    if (i % 100 == 0) {
      double phi = cal_order_para(birds);
      auto t_now = std::chrono::system_clock::now();
      std::chrono::duration<double> elapsed_seconds = t_now - t_start;
      fout_phi << i << "\t" << phi << "\t" << elapsed_seconds.count() << "s\n";
    }
    if (flag_cg && ((i < 20) ||
      (i < 80 && i % 2 == 0) ||
      (i < 160 && i % 4 == 0) ||
      (i < 400 && i % 8 == 0) ||
      (i < 1000 && i % 16 == 0) ||
      (i % 32 == 0))) {
      output_coarse_grained_snap(i, n_cg, vx_cg, vy_cg, n_boxes,
        box_size_cg, fout_snap, pos_idx_pair, birds);
    }
  }
}

void end() {
  auto t_end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = t_end - t_start;
  std::time_t end_time = std::chrono::system_clock::to_time_t(t_end);
  std::cout << "Finished computation at " << std::ctime(&end_time)
    << "elapsed time: " << elapsed_seconds.count() << "s\n";
}
