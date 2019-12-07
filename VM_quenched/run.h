#pragma once
#include <vector>
#include "vect.h"
#include "comn.h"
#include "output.h"
#include "configure.h"
#ifdef USE_MPI
#include "mpi.h"
#endif

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
  Vec_2<double> sum(0., 0.);
  for (auto &i: rand_field) {
    sum += i;
  }
  std::cout << "sum of the random fields: " << sum << std::endl;
}

template <typename TNode, typename TCellList, typename TPairForce, typename TMove>
void run(double eta, double eps, const Vec_2<double> &gl_l, unsigned long long seed, int n_step,
         TPairForce pair_force, TMove move,
         std::vector<TNode> &p_arr,
         TCellList &cl) {
  ini_order_para_exporter(eta, eps, gl_l.x, seed);
#ifndef USE_MPI
  for (int i = 1; i <= n_step; i++) {
    cl.for_each_pair(pair_force);
    move();
    if (i % 100 == 0) {
      double phi, theta;
      cal_order_para(p_arr, phi, theta);
      output_order_para(phi, theta);
      if (i % 100000 == 0) {
        std::cout << "t = " << i << std::endl;
        output_coarse_grained_snap(p_arr, i, gl_l, eta, eps, seed);
      }
    }
  }
#else
  int i_step = 1;
  MPI_Win win;
  int finished_proc = 0;
  int proc_size;
  MPI_Win_create(&i_step, 1, sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
  MPI_Comm_size(MPI_COMM_WORLD, &proc_size);
  while (finished_proc < proc_size) {
    cl.for_each_pair(pair_force);
    move();
    if (i_step % 100 == 0) {
      double phi, theta;
      cal_order_para(p_arr, phi, theta);
      output_order_para(phi, theta);
      if (i_step % 100000 == 0) {
        std::cout << "t = " << i_step << std::endl;
        output_coarse_grained_snap(p_arr, i_step, gl_l, eta, eps, seed);
      }
    }
    i_step++;
    if (i_step > n_step && i_step % 1000 == 1) {
      finished_proc = 0;
      for (int j = 0; j < proc_size; j++) {
        int remote_i_step;
        MPI_Win_lock(MPI_LOCK_SHARED, j, 0, win);
        MPI_Get(&remote_i_step, 1, MPI_INT, j, 0, 1, MPI_INT, win);
        MPI_Win_unlock(j, win);
        if (remote_i_step > n_step) {
          finished_proc++;
        }
      }
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}