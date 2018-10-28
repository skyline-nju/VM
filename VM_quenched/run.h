#pragma once
#include <vector>
#include "vect.h"
#include "comn.h"
#include "output.h"
#include "configure.h"

void set_para(int argc, char* argv[]);

template <typename TPar, typename TInteg>
void integrate_all(const std::vector<TPar> &p_arr, TInteg my_integ) {
  auto end = p_arr.end();
  for (auto it = p_arr.begin(); it != p_arr.end(); ++it) {
    my_integ(*it);
  }
}

template <typename TRan>
void ini_rand_torque(std::vector<double> &rand_torque, double eps, int n, TRan &myran) {
  rand_torque.reserve(n);
  double d = 1.0 / (n - 1);
  for (int i = 0; i < n; i++) {
    rand_torque.push_back((-0.5 + i * d) * eps * 2.0 * PI);
  }
  shuffle(rand_torque, myran);
  double sum = 0;
  for (auto i : rand_torque) {
    sum += i;
  }
  std::cout << "sum of the random torques: " << sum << std::endl;
}

template <typename TRan>
void ini_rand_field(std::vector<Vec_2<double>> &rand_field, double eps, int n, TRan &myran) {
  const double dtheta = PI * 2. / n;
  double *theta_arr = new double[n];
  for (int i = 0; i < n; i++) {
    theta_arr[i] = i * dtheta;
  }
  shuffle(theta_arr, n, myran);
  rand_field.reserve(n);
  for (int i = 0; i < n; i++) {
    rand_field.emplace_back(eps * cos(theta_arr[i]), eps * sin(theta_arr[i]));
  }
  delete[] theta_arr;
}

template <typename TNode, typename TCellList, typename TPairForce, typename TMove>
void run(double eta, double eps, const Vec_2<double> &gl_l, int seed, int n_step,
         TPairForce pair_force, TMove move,
         std::vector<TNode> &p_arr,
         TCellList &cl) {
  //ini_order_para_exporter(eta, eps, gl_l.x, seed);
  double phi_mean = 0;
  int count = 0;
  for (int i = 1; i <= n_step; i++) {
    cl.for_each_pair(pair_force);
    move();
    if (i % 100 == 0) {
      double phi, theta;
      cal_order_para(p_arr, phi, theta);
      //output_order_para(phi, theta);
      if (i >= 100000) {
        phi_mean += phi;
        count++;
      }
      if (i % 10000 == 0) {
        std::cout << "t = " << i << std::endl;
        //if (i % 100000 == 0) {
        //  output_coarse_grained_snap(p_arr, i, gl_l, eta, eps, seed);
        //}
      }
    }
  }
  std::cout << "time-averaged order parameter: " << phi_mean / count << std::endl;
}

template <typename TNode, typename TCellList, typename TPairForce, typename TMove>
void run_MSD(double eta, double eps, const Vec_2<double> &gl_l, int seed, int n_step,
             TPairForce pair_force, TMove move, std::vector<TNode> &p_arr, TCellList &cl) {
  std::vector<int> t0_arr;
  int t0_sep = 500000;
  int t_eq = 100000;
  int t = 100000;
  int record_sep = 100;
  while (t < n_step / 2) {
    t0_arr.push_back(t);
    t += t0_sep;
  }
  std::vector<MSD> msd_arr;
  msd_arr.reserve(t0_arr.size());
  ini_order_para_exporter(eta, eps, gl_l.x, seed);
  ini_XY_exporter(eta, eps, gl_l.x, seed);
  ini_traj(eta, eps, gl_l.x, seed);
  for (int i = 1; i <= n_step; i++) {
    cl.for_each_pair(pair_force);
    move();
    if (i % 100 == 0) {
      double phi, theta;
      cal_order_para(p_arr, phi, theta);
      output_order_para(phi, theta);
      if (i > t_eq) {
        output_traj(p_arr, i, gl_l);
      }
    }
    for (auto t0: t0_arr) {
      if (i == t0) {
        msd_arr.emplace_back(p_arr, t0, gl_l, record_sep);
      }
    }
    for (auto &msd: msd_arr) {
      msd.record(p_arr, i, gl_l);
    }
    if (i >= t_eq && i <= t_eq + 100) {
      output_XY(p_arr, i, gl_l);
    }
    if (i % t_eq == 0) {
      std::cout << "t = " << i << std::endl;
    }
  }
  output_mean_msd(msd_arr, eta, eps, gl_l.x, seed);
}