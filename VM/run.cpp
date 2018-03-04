#include "run.h"
using namespace std;

void ini_birds(Node **p_bird, Ran *myran, const cmdline::parser &cmd) {
  if (!cmd.exist("rho0") && !cmd.exist("file")) {
    cerr << cmd.usage();
    exit(1);
  } else if (cmd.exist("file")) {
    InSnapshot isnap(cmd.get<string>("file"));
    vector<float> x;
    vector<float> y;
    vector<float> theta;
    double Lx;
    double Ly;
    isnap.from_file(cmd.get<int>("idx_frame"), Lx, Ly, x, y, theta);
    Node::rho_0 = x.size() / (Lx*Ly);
    *p_bird = Node::ini_from_snap(Lx, Ly, x, y, theta);
    cout << "rho_0 = " << Node::rho_0 << endl;
  } else {
    Node::rho_0 = cmd.get<double>("rho0");
    if (cmd.get<string>("ini_mode") == "left") {
      *p_bird = Node::ini_move_left(myran);
    } else if (cmd.get<string>("ini_mode") == "rand") {
      *p_bird = Node::ini_rand(myran);
    } else {
      cerr << cmd.usage();
      exit(1);
    }
  }
}

void ini_rand_torques(double **disorder, int n, double epsilon,
                      unsigned long long seed) {
  if (epsilon > 0) {
    Ran myran(seed);
    double d = 1.0 / (n - 1);
    double *torque = new double[n];
    for (int i = 0; i < n; i++)
      torque[i] = (-0.5 + i *d) * epsilon * 2.0 * PI;
    shuffle(torque, n, &myran);
    *disorder = torque;
  }
}


void update_coor(Node *bird, Ran* myran, double eta, const double *disorder, bool vicskeShake) {
  static double eta2PI = eta * 2 * PI;

  //calculating noise
  double *noise = new double[Node::N];
  if (disorder) {
    for (int i = 0; i < Node::N; i++)
      noise[i] = (myran->doub() - 0.5) * eta2PI + disorder[bird[i].cell_idx];
  } else {
    for (int i = 0; i < Node::N; i++)
      noise[i] = (myran->doub() - 0.5) * eta2PI;
  }

  //updating coordination
  if (vicskeShake) {
    for (int i = 0; i < Node::N; i++) {
      bird[i].move(noise[i], myran);
    }
  } else {
    for (int i = 0; i < Node::N; i++) {
      bird[i].move(noise[i]);
    }

  }
  delete[] noise;
  noise = nullptr;
}

void run(Node *bird, Grid *cell, Ran *myran, int nStep,
	       double eta, const double *disorder, Output &out, bool vicsekShake) {
  for (int i = 1; i <= nStep; i++) {
    Grid::all_pairs(cell);
    update_coor(bird, myran, eta, disorder, vicsekShake);
    Grid::refresh(cell, bird);
    out.out(bird, Node::N, i);
  }
}

void run_raw(Node * bird, Grid * cell, Ran * myran, int nStep,
						 double eta, const double * disorder, bool vicsekShake) {
	for (int i = 1; i <= nStep; i++) {
		Grid::all_pairs(cell);
		update_coor(bird, myran, eta, disorder, vicsekShake);
		Grid::refresh(cell, bird);
	}
}
