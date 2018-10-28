#include "output.h"
#include <fstream>
#include <iomanip>

std::ofstream fout_phi;
std::ofstream fout_XY;
std::ofstream fout_traj;

#ifdef RAND_FIELD
char disorder_t[] = "rf";
#else
char disorder_t[] = "rt";
#endif

void ini_order_para_exporter(double eta, double eps, double L, unsigned long long seed) {
  char filename[100];
  snprintf(filename, 100, "phi_%s_%g_%g_%g_%llu.dat", disorder_t, L, eta, eps, seed);
  fout_phi.open(filename);
}

void output_order_para(double phi, double theta) {
  fout_phi << std::setprecision(8) << std::fixed << phi << "\t" << theta << std::endl;
}

void ini_XY_exporter(double eta, double eps, double L, unsigned long long seed) {
  char filename[100];
  snprintf(filename, 100, "%s_%g_%g_%g_%llu.extxyz", disorder_t, L, eta, eps, seed);
  fout_XY.open(filename);
}

void output_XY_frame_head(const Vec_2<double>& l, int t, int n_par) {
  fout_XY << n_par << "\n";
  fout_XY << "Lattice=\"" << l.x << " 0 0 0 " << l.y << " 0 0 0 1\" "
          << "Properties=species:S:1:pos:R:2:mass:M:1 "
          << "Time=" << t;
}

void output_XY_frame_tail() {
  fout_XY << std::endl;
}

void output_XY_one(double x, double y, double theta) {
  fout_XY << "\n" << std::fixed << std::setprecision(3) << "N\t"
    << x << "\t" << y << "\t" << theta;
}

void ini_traj(double eta, double eps, double L, unsigned long long seed) {
  char filename[100];
  snprintf(filename, 100, "traj_%s_%g_%g_%g_%llu.dat", disorder_t, L, eta, eps, seed);
  fout_traj.open(filename);
}

void output_traj(const std::vector<double>& x_arr, const std::vector<double>& y_arr, int t) {
  fout_traj << std::fixed << std::setprecision(3) << t;
  int n = x_arr.size();
  for (int i = 0; i < n; i++) {
    fout_traj << "\t" << x_arr[i] << "\t" << y_arr[i];
  }
  fout_traj << std::endl;
}

void output_coarse_grained_snap(unsigned short* num, float* vx, float* vy, int n, int t,
                                double eta, double eps, double L, unsigned long long seed) {
  char filename[100];
  snprintf(filename, 100, "snap/cg%s_%g_%g_%g_%llu_%07d.bin", disorder_t, L, eta, eps, seed, t);
  std::ofstream fout(filename, std::ios::binary);
  fout.write((char*)&num[0], sizeof(unsigned short) * n);
  fout.write((char*)&vx[0], sizeof(float) * n);
  fout.write((char*)&vy[0], sizeof(float) * n);
  fout.close();
}

void MSD::output(double eta, double eps, double L, unsigned long long seed) {
  char filename[100];
  snprintf(filename, 100, "msd_%s_%g_%g_%g_%llu.dat", disorder_t, L, eta, eps, seed);
  std::ofstream fout(filename);
  int n = t_arr_.size();
  for (int i = 0; i < n; i++) {
    fout << t_arr_[i] << "\t" << r_square_[i] << "\n";
  }
  fout.close();
}

void output_mean_msd(const std::vector<MSD>& msd_arr,
                     double eta, double eps, double L, unsigned long long seed) {
  int n_frame = msd_arr[0].r_square().size();
  int n_msd = msd_arr.size();
  std::vector<double> mean_msd(msd_arr[0].r_square());
  std::vector<int> count_arr(n_frame, 1);
  for (int i = 1; i < n_msd; i++) {
    int len = msd_arr[i].r_square().size();
    for (int j = 0; j < len; j++) {
      mean_msd[j] += msd_arr[i].r_square()[j];
      count_arr[j]++;
    }
  }
  for (int i = 0; i < n_frame; i++) {
    mean_msd[i] /= count_arr[i];
  }

  char filename[100];
#ifdef RAND_FIELD
  char disorder_t[] = "rf";
#else
  char disorder_t[] = "rt";
#endif
  snprintf(filename, 100, "msd_%s_%g_%g_%g_%llu.dat", disorder_t, L, eta, eps, seed);
  std::ofstream fout(filename);
  for (int i = 0; i < n_frame; i++) {
    fout << msd_arr[0].t_arr()[i] << "\t" << mean_msd[i] << "\n";
  }
  fout.close();
}
