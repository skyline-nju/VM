#ifndef IO_DATA_H
#define IO_DATA_H
#include <fstream>
#include <ctime>
#include <chrono>
#include "comn.h"
#include "cmdline.h"
#include "vm.h"
#include "gsd.h"

void ini_output(const cmdline::parser &cmd, const VM *birds, const std::string& prefix);

void output(int i, const VM* birds);

class _Writer {
public:
  _Writer(const cmdline::parser &cmd, const std::string& prefix) : idx_frame(0) {}
  virtual void set_frames(const cmdline::parser &cmd) = 0;
  virtual void write(int i, const VM *bird) = 0;

protected:
  std::vector<int> frames;
  int idx_frame;
};

class OrderParaWriter : public _Writer {
public:
  OrderParaWriter(const cmdline::parser &cmd, const std::string& prefix);
  void set_frames(const cmdline::parser &cmd);
  void write(int i, const VM* bird);

private:
  std::ofstream fout_; 
};

class LogWriter : public _Writer {
public:
  LogWriter(const cmdline::parser &cmd, const std::string& prefix);
  void set_frames(const cmdline::parser &cmd);
  void write(int i, const VM *bird);
private:
  std::chrono::time_point<std::chrono::system_clock> t_start;
  std::ofstream fout_;
};

class CoarseGrainSnapWriter : public _Writer {
public:
  CoarseGrainSnapWriter(const cmdline::parser &cmd, const std::string& prefix);
  void set_frames(const cmdline::parser &cmd);
  void write(int i, const VM* bird);

private:
  double l;
  int ncols;
  int ncells;
  double dt;
  bool is_continous;
  std::ofstream fout_;
};

template <typename T1, typename T2>
void coarse_grain(T1 *num_cg, T2 *vx_cg, T2 *vy_cg, double &vx_m, double &vy_m,
                  double l, int ncols, int ntot, const VM *birds) {
  for (int i = 0; i < ntot; i++) {
    num_cg[i] = 0;
    vx_cg[i] = 0;
    vy_cg[i] = 0;
  }
  vx_m = 0;
  vy_m = 0;
  int N = birds->get_num_birds();
  for (int i = 0; i < N; i++) {
    double x, y, vx, vy;
    birds->get_x(i, x, y);
    birds->get_v(i, vx, vy);
    int col = int(x / l);
    int row = int(y / l);
    int j = col + ncols * row;
    num_cg[j]++;
    vx_cg[j] += vx;
    vy_cg[j] += vy;
    vx_m += vx;
    vy_m += vy;
  }
  vx_m /= N;
  vy_m /= N;
}



/*******************************************************************************************
 * @brief Basic class for exporting data.
 *
 * Define the timming to dump data.
 ******************************************************************************************/
class ExporterBase {
public:
  ExporterBase() : n_step_(0) {}

  ExporterBase(int start, int n_step, int sep) : start_(start), n_step_(n_step), sep_(sep) {
    set_lin_frame(start, n_step, sep);
  }

  void set_lin_frame(int start, int n_step, int sep);

  bool need_export(const int i_step);

protected:
  int n_step_;    // total steps to run
  int sep_;
  int start_ = 0; // The first step 
private:
  std::vector<int> frames_arr_; // frames that need to export
  std::vector<int>::iterator frame_iter_;
};

/*******************************************************************************************
 * @brief Read and output snapshot in gsd format
 ******************************************************************************************/
class Snap_GSD_2 : public ExporterBase {
public:
  Snap_GSD_2(const std::string& filename, int start, int n_step, int sep,
             double Lx, double Ly, const std::string& open_flag);

  ~Snap_GSD_2();

  template <typename TPar>
  void get_data_from_par(const std::vector<TPar>& p_arr, float* pos);

  void get_data_from_par(const VM* birds,
                         float* pos, uint32_t* type_id=nullptr) const;

  uint64_t get_time_step();

  template <typename TPar>
  void dump(int i_step, const std::vector<TPar>& p_arr);
  void dump(int i_step, const VM*birds);

  template <typename TPar>
  void read(int i_frame, std::vector<TPar>& p_arr);
  //TODO

  template <typename TPar>
  void read_last_frame(std::vector<TPar>& p_arr);
  //TODO

private:
  gsd_handle* handle_ = nullptr;
  double half_Lx_;
  double half_Ly_;
};

template <typename TPar>
void Snap_GSD_2::get_data_from_par(const std::vector<TPar>& p_arr, float* pos
  //, uint32_t* type_id
) {
  size_t n_par = p_arr.size();
  for (size_t j = 0; j < n_par; j++) {
    size_t j3 = j * 3;
    pos[j3    ] = p_arr[j].x - half_Lx_;
    pos[j3 + 1] = p_arr[j].y - half_Ly_;
    pos[j3 + 2] = p_arr[j].get_theta();
    //type_id[j] = p_arr[j].type_id;
  }
}


template<typename TPar>
void Snap_GSD_2::dump(int i_step, const std::vector<TPar>& p_arr) {
  if (need_export(i_step)) {
    uint32_t n_par = p_arr.size();
    float* pos = new float[n_par * 3];
    //uint32_t *type_id = new uint32_t[n_par];
    get_data_from_par(p_arr, pos);
    uint64_t step = get_time_step();
  
    std::cout << "dump frame " << step << std::endl;
    gsd_write_chunk(handle_, "configuration/step", GSD_TYPE_UINT64, 1, 1, 0, &step);
    gsd_write_chunk(handle_, "particles/N", GSD_TYPE_UINT32, 1, 1, 0, &n_par);
    gsd_write_chunk(handle_, "particles/position", GSD_TYPE_FLOAT, n_par, 3, 0, pos);
    //gsd_write_chunk(handle_, "particles/typeid", GSD_TYPE_UINT32, n_par, 1, 0, type_id);
    gsd_end_frame(handle_);
    delete[] pos;
    //delete []type_id;
  }
}


template <typename TPar>
void Snap_GSD_2::read(int i_frame, std::vector<TPar>& p_arr) {
  uint32_t n_par;
  const gsd_index_entry* chunk = gsd_find_chunk(handle_, i_frame, "particles/N");
  gsd_read_chunk(handle_, &n_par, chunk);
  std::cout << "frame " << i_frame  <<": find " << n_par << " particles" << std::endl;

  float* pos = new float[n_par * 3];
  chunk = gsd_find_chunk(handle_, i_frame, "particles/position");
  gsd_read_chunk(handle_, pos, chunk);
  //uint32_t* type_id = new uint32_t[n_par];
  //chunk = gsd_find_chunk(handle_, i_frame, "particles/typeid");
  //gsd_read_chunk(handle_, type_id, chunk);

  p_arr.reserve(n_par);

  double Lx = half_Lx_ * 2;
  double Ly = half_Ly_ * 2;

  for (int j = 0; j < n_par; j++) {
    size_t j3 = j * 3;
    double x = pos[j3] + half_Lx_;
    double y = pos[j3 + 1] + half_Ly_;
    double theta = pos[j3 + 2];
    
    if (x < 0) {
      x += Lx;
    } else if (x >= Lx) {
      x -= Lx;
    }
    if (y < 0) {
      y += Ly;
    } else if (y >= Ly) {
      y -= Ly;
    }

    p_arr.emplace_back(x, y, theta);
  }

  delete[] pos;
  //delete[] type_id;
}


template <typename TPar>
void Snap_GSD_2::read_last_frame(std::vector<TPar>& p_arr) {
  int nframes = gsd_get_nframes(handle_);
  if (nframes < 1) {
    std::cout << "Error, nframes=" << nframes << std::endl;
    exit(1);
  } else {
    read(nframes-1, p_arr);
  }
}


#endif
