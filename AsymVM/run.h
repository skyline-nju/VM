#ifndef RUN_H
#define RUN_H
#include <cstring>
#include <iostream>

#include "grid.h"
//#include "exporter.h"


template<typename TRan, typename TPar>
void run_one_step(std::vector<TPar>& p_arr, TRan& myran, Grid<TPar>& cells,
                  double eta, double Lx, double Ly, double v0) {

  cells.all_pairs();
  for (auto& p : p_arr) {
    p.move(eta, myran, Lx, Ly, v0);
  }
  cells.refresh(p_arr);
}


template<typename TRan, typename TPar>
void ini_rand(std::vector<TPar>& p_arr, int N, double Lx, double Ly, TRan& myran) {
  for (int i = 0; i < N; i++) {
    p_arr.emplace_back(myran, Lx, Ly);
  }
}

template<typename TRan, typename TPar, typename TSnap>
void ini_particles(std::vector<TPar>& p_arr, TRan& myran, Grid<TPar>& cells,
                   int n_par, double Lx, double Ly, const std::string& ini_mode, TSnap& snap) {
  p_arr.reserve(n_par);
  if (ini_mode == "rand") {
    ini_rand(p_arr, n_par, Lx, Ly, myran);
  } else if (ini_mode == "resume") {
    snap.read_last_frame(p_arr);
  }
  cells.link_nodes(p_arr);
}

#endif
