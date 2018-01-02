#include "vm.h"

VM::VM(const cmdline::parser &cmd, Ran &myran) {
  Lx = cmd.get<double>("Lx");
  Ly = cmd.exist("Ly") ? cmd.get<double>("Ly") : Lx;
  v0 = cmd.get<double>("v0");
  double rho0 = cmd.get<double>("rho0");
  N = int(Lx * Ly * rho0);
  eta = cmd.get<double>("eta");
  extra_torque = cmd.get<double>("ext_torque") * PI;
  double eps = cmd.get<double>("eps");
  if (eps > 0) {
    int tot = int(Lx * Ly);
    random_torque = new double[tot];
    double d = 1.0 / (tot - 1);
    for (int i = 0; i < tot; i++)
      random_torque[i] = (-0.5 + i * d) * eps * 2 * PI;
    myran.shuffle(random_torque, tot);
  } else {
    random_torque = nullptr;
  }
  if (cmd.exist("dt")) {
    double dt = cmd.get<double>("dt");
    V_conti::h = dt;
    V_conti::sqrt_24_Dr_h = std::sqrt(24 * eta * dt);
  }

#ifdef _MSC_VER
  std::cout << "System size: " << Lx << " x " << Ly << "\n";
  std::cout << "eta = " << eta << "\n";
  std::cout << "eps = " << eps << "\n";
  std::cout << "extra torque = " << extra_torque << "\n";
#endif
}

VM::~VM() {
  if (random_torque) {
    delete[] random_torque;
    random_torque = nullptr;
  }
}

void VM::output_data(double * x, double * y, double * vx, double * vy) const {
  for (int i = 0; i < N; i++) {
    get_x(i, x[0], y[0]);
    get_v(i, vx[0], vy[0]);
  }
}

void VM::ini(const cmdline::parser & cmd, Ran & myran) {
  if (cmd.exist("defect_sep")) {
    double d = cmd.get<double>("defect_sep");
    double alpha = cmd.get<int>("defect_mode") * PI * 0.5;
    create_defect_pair(myran, d, alpha);
#ifdef _MSC_VER
    std::cout << "create a pair of defects with separation of "
      << d * 2 << "\n";
#endif
  } else {
    create_random(myran);
#ifdef _MSC_VER
    std::cout << "create particles with random positions and velocities\n";
#endif
  }
}

void VM::create_random(Ran & myran) {
  double *x = new double[N];
  double *y = new double[N];
  double *vx = new double[N];
  double *vy = new double[N];
  for (int i = 0; i < N; i++) {
    x[i] = myran.doub() * Lx;
    y[i] = myran.doub() * Ly;
    myran.circle_point_picking(vx[i], vy[i]);
  }
  input_data(x, y, vx, vy);
  delete[] x;
  delete[] y;
  delete[] vx;
  delete[] vy;
}

void VM::create_defect_pair(Ran &myran,
                            double defect_sep,
                            double rot_angle) {
  double *x = new double[N];
  double *y = new double[N];
  double *vx = new double[N];
  double *vy = new double[N];
  for (int i = 0; i < N; i++) {
    x[i] = myran.doub() * Lx;
    y[i] = myran.doub() * Ly;
    double X = x[i] - 0.5 * Lx;
    double Y = y[i] - 0.5 * Ly;
    double phi1 = std::atan2(Y, X - defect_sep);
    double phi2 = std::atan2(Y, X + defect_sep);
    double phi = phi1 - phi2 + rot_angle;
    vx[i] = std::cos(phi);
    vy[i] = std::sin(phi);
  }
  input_data(x, y, vx, vy);
  delete[] x;
  delete[] y;
  delete[] vx;
  delete[] vy;
}

void VM::get_v_mean(double & vx_m, double & vy_m) const {
  vx_m = 0;
  vy_m = 0;
  for (int i = 0; i < N; i++) {
    double vx, vy;
    get_v(i, vx, vy);
    vx_m += vx;
    vy_m += vy;
  }
  vx_m /= N;
  vy_m /= N;
}