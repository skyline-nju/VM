#include "exporter.h"
#include <netcdf.h>

double rho0;
double eta;
double eps;
unsigned long long seed;
double Lx;
double Ly;
int n_par;
int n_step;
std::string folder;
std::string base_name;
int deflate_level = 6;

void set_output(const cmdline::parser& cmd,
                LogExporter** log_ex,
                OrderParaExpoter** order_ex,
                CoarseGrainSnapExporter** cg_ex,
                SpatialCorr2dExporter** corr_ex) {
  if (cmd.exist("output")) {
    Lx = cmd.get<double>("Lx");
    Ly = cmd.exist("Ly") ? cmd.get<double>("Ly") : Lx;

    rho0 = cmd.get<double>("rho0");
    eta = cmd.get<double>("eta");
    eps = cmd.get<double>("eps");
    n_par = Lx * Ly * rho0;
    n_step = cmd.get<int>("nstep");
    seed = cmd.get<unsigned long long>("seed");

    folder = cmd.get<std::string>("output");
    mkdir(folder);
    folder += delimiter;
    char str[100];
    snprintf(str, 100, "%g_%g_%g_%g_%g_%llu", eta, eps, rho0, Lx, Ly, seed);
    base_name = str;

    if (cmd.exist("log_dt")) {
      *log_ex = new LogExporter(cmd);
    } else {
      *log_ex = nullptr;
    }
    if (cmd.exist("order_dt")) {
      *order_ex = new OrderParaExpoter(folder + "p_" + base_name + ".dat",
                                       n_step, cmd.get<int>("order_dt"));
    } else {
      *order_ex = nullptr;
    }

    int ncols = cmd.exist("cg_ncols") ? cmd.get<int>("cg_ncols") : int(Lx);
    int nrows = cmd.exist("cg_nrows") ? cmd.get<int>("cg_nrows") : int(Ly);

    if (cmd.exist("cg_dt")) {
      int end = cmd.exist("cg_beg") ? cmd.get<int>("cg_end") : n_step;
      *cg_ex = new CoarseGrainSnapExporter(folder + "cg_" + base_name + ".nc",
                                          ncols, nrows, end,
                                          cmd.get<int>("cg_dt"),
                                          cmd.get<int>("cg_beg"));
    } else {
      *cg_ex = nullptr;
    }
    if (cmd.exist("corr_dt")) {
      *corr_ex = new SpatialCorr2dExporter(folder + "cr_" + base_name + ".nc",
                                           ncols, nrows, n_step,
                                           cmd.get<int>("corr_dt"),
                                           cmd.get<int>("corr_beg"));
    } else {
      *corr_ex = nullptr;
    }
  }
}

/*************************************************************************//**
 * @brief Construct a new Log Exporter:: Log Exporter object
 * 
 * @param cmd Cmdline parser
 ****************************************************************************/
LogExporter::LogExporter(const cmdline::parser& cmd)
  : BaseLogExporter(folder + base_name + ".log", n_par,
                    n_step, cmd.get<int>("log_dt")) {
  fout_ << "\n-------- Parameters --------";
  fout_ << "\nParticle number=" << n_par;
  fout_ << "\nrho0=" << rho0;
  fout_ << "\nLx=" << Lx;
  fout_ << "\nLy=" << Ly;
  fout_ << "\neta=" << eta;
  fout_ << "\nepsilon=" << eps;
  fout_ << "\nn_step=" << n_step;;
  fout_ << "\nseed=" << seed;
  fout_ << "\n\n-------- RUN --------";
  fout_ << "\ntime step\telapsed time" << std::endl;
}

OrderParaExpoter::OrderParaExpoter(const std::string & filename,
                                   int end, int sep, int start)
  :BaseExporter(end, sep, start), fout_(filename) {
  std::cout << "begin to record the order parameter\n";
}

// check whether there is error when outputting netcdf file
void check_err(const int stat, const int line, const char * file) {
  if (stat != NC_NOERR) {
    (void)fprintf(stderr, "line %d of %s: %s\n", line, file, nc_strerror(stat));
    fflush(stderr);
    exit(1);
  }
}

/**
* \brief output global parameters into netcdf file
* \param ncid  Id for netcdf file
*/
void put_global_para(int ncid) {
  auto stat = nc_put_att_text(ncid, NC_GLOBAL, "title", 21, "quenched Vicsek model");
  check_err(stat, __LINE__, __FILE__);
  stat = nc_put_att_double(ncid, NC_GLOBAL, "density", NC_DOUBLE, 1, &rho0);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_put_att_double(ncid, NC_GLOBAL, "Lx", NC_DOUBLE, 1, &Lx);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_put_att_double(ncid, NC_GLOBAL, "Ly", NC_DOUBLE, 1, &Ly);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_put_att_double(ncid, NC_GLOBAL, "eta", NC_DOUBLE, 1, &eta);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_put_att_double(ncid, NC_GLOBAL, "epsilon", NC_DOUBLE, 1, &eps);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_put_att_ulonglong(ncid, NC_GLOBAL, "seed", NC_UINT64, 1, &seed);
  check_err(stat, __LINE__, __FILE__);
}

/*************************************************************************//**
 * @brief Construct a new Coarse Grain Snap Exporter:: Coarse Grain Snap Exporter object
 * 
 * @param filename  File to output
 * @param ncols     number of cols
 * @param nrows     number of rows
 * @param end       Last time step to record
 * @param sep       Spacing between two frames
 * @param start     First time step to record
 ****************************************************************************/
CoarseGrainSnapExporter::CoarseGrainSnapExporter(const std::string & filename, // NOLINT
                                                 int ncols, int nrows,
                                                 int end, int sep, int start)
  :BaseExporter(end, sep, start), frame_len_(NC_UNLIMITED), time_idx_{0},
  ncols_(ncols), nrows_(nrows), lx_(Lx / ncols), ly_(Ly / nrows) {
  std::cout << lx_ << "\t" << ly_ << "\n";
  auto stat = nc_create(filename.c_str(), NC_NETCDF4, &ncid_);
  check_err(stat, __LINE__, __FILE__);

  /* dimension ids */
  int frame_dim;
  int row_dim;
  int col_dim;

  /* define dimensions */
  stat = nc_def_dim(ncid_, "frame", frame_len_, &frame_dim);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid_, "nrows", nrows_, &row_dim);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid_, "ncols", ncols_, &col_dim);
  check_err(stat, __LINE__, __FILE__);

  /* define variables */
  int time_dims[1] = {frame_dim};
  stat = nc_def_var(ncid_, "time", NC_INT, 1, time_dims, &time_id_);
  check_err(stat, __LINE__, __FILE__);
  int corr_dims[3] = {frame_dim, row_dim, col_dim};
  stat = nc_def_var(ncid_, "num", NC_SHORT, 3, corr_dims, &num_id_);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_def_var(ncid_, "vx", NC_FLOAT, 3, corr_dims, &vx_id_);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_def_var(ncid_, "vy", NC_FLOAT, 3, corr_dims, &vy_id_);
  check_err(stat, __LINE__, __FILE__);

  /* set chunk and deflating */
  {
    // time and mean_velocity
    size_t time_block = get_n_frames() / 10;
    if (time_block >= 4096)
      time_block = 4096;
    size_t chunk[1] = {time_block};
    stat = nc_def_var_chunking(ncid_, time_id_, NC_CHUNKED, chunk);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_def_var_deflate(ncid_, time_id_, NC_SHUFFLE, 1, deflate_level);
    check_err(stat, __LINE__, __FILE__);
  }
  {
    // coarse-grained snapshots
    size_t time_block = 4096 * 12 * 2 / (ncols * nrows);
    if (time_block < 1)
      time_block = 1;
    else if (time_block > get_n_frames())
      time_block = get_n_frames();
    size_t chunk[3] = {time_block, nrows_, ncols_};
    stat = nc_def_var_chunking(ncid_, num_id_, NC_CHUNKED, chunk);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_def_var_deflate(ncid_, num_id_, NC_SHUFFLE, 1, deflate_level);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_def_var_chunking(ncid_, vx_id_, NC_CHUNKED, chunk);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_def_var_deflate(ncid_, vx_id_, 0, 1, deflate_level);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_def_var_chunking(ncid_, vy_id_, NC_CHUNKED, chunk);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_def_var_deflate(ncid_, vy_id_, 0, 1, deflate_level);
    check_err(stat, __LINE__, __FILE__);
  }

  /* put global values */
  put_global_para(ncid_);

  /* leave define mode */
  stat = nc_enddef(ncid_);
  check_err(stat, __LINE__, __FILE__);
}

void CoarseGrainSnapExporter::write_frame(int i_step, const short* num,
                                          const float* vx, const float* vy) {
  /* time step */
  {
    auto stat = nc_put_var1(ncid_, time_id_, time_idx_, &i_step);
    check_err(stat, __LINE__, __FILE__);
  }
  /* coarse-grained snapshots */
  {
    size_t startset[3] = {time_idx_[0], 0, 0};
    size_t countset[3] = {1, nrows_, ncols_};
    auto stat = nc_put_vara(ncid_, num_id_, startset, countset, &num[0]);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_put_vara(ncid_, vx_id_, startset, countset, &vx[0]);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_put_vara(ncid_, vy_id_, startset, countset, &vy[0]);
    check_err(stat, __LINE__, __FILE__);
  }
  ++time_idx_[0]; //! update time frame
}

/*************************************************************************//**
 * @brief Construct a new Spatial Corr 2d Exporter:: Spatial Corr 2d Exporter object
 * 
 * @param filename  Filename to output
 * @param ncols     number of cols
 * @param nrows     number of rows
 * @param end       Last time step to record
 * @param sep       Spacing between two frames
 * @param start     First time step to record
 ****************************************************************************/
SpatialCorr2dExporter::SpatialCorr2dExporter(const std::string& filename, // NOLINT
                                             int ncols, int nrows,
                                             int end, int sep, int start)
  : SpatialCorr2d(ncols, nrows, Lx, Ly), BaseExporter(end, sep, start),
  frame_len_(NC_UNLIMITED), time_idx_{0} {

  auto stat = nc_create(filename.c_str(), NC_NETCDF4, &ncid_);
  check_err(stat, __LINE__, __FILE__);

  /* dimension ids */
  int frame_dim;
  int row_dim;
  int col_dim;
  int spatial_dim;

  /* define dimensions */
  stat = nc_def_dim(ncid_, "frame", frame_len_, &frame_dim);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid_, "nrows", nrows_ / 2, &row_dim); //! only half rows are needed
  check_err(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid_, "ncols", ncols_, &col_dim);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid_, "spatial", 2, &spatial_dim);
  check_err(stat, __LINE__, __FILE__);

  /* define variables */
  int time_dims[1] = {frame_dim};
  stat = nc_def_var(ncid_, "time", NC_INT, 1, time_dims, &time_id_);
  check_err(stat, __LINE__, __FILE__);
  int corr_dims[3] = {frame_dim, row_dim, col_dim};
  stat = nc_def_var(ncid_, "C_rho", NC_FLOAT, 3, corr_dims, &corr_rho_id_);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_def_var(ncid_, "C_v", NC_FLOAT, 3, corr_dims, &corr_v_id_);
  check_err(stat, __LINE__, __FILE__);
  int v_mean_dims[2] = {frame_dim, spatial_dim};
  stat = nc_def_var(ncid_, "mean_velocity", NC_DOUBLE, 2, v_mean_dims, &v_mean_id_);
  check_err(stat, __LINE__, __FILE__);

  /* set chunk and deflating */
  {
    // time and mean_velocity
    size_t time_block = get_n_frames() / 10;
    if (time_block >= 4096)
      time_block = 4096;
    size_t chunk[1] = {time_block};
    stat = nc_def_var_chunking(ncid_, time_id_, NC_CHUNKED, chunk);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_def_var_deflate(ncid_, time_id_, NC_SHUFFLE, 1, deflate_level);
    check_err(stat, __LINE__, __FILE__);
    size_t chunk2[2] = {time_block, 2};
    stat = nc_def_var_chunking(ncid_, v_mean_id_, NC_CHUNKED, chunk2);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_def_var_deflate(ncid_, v_mean_id_, 0, 1, deflate_level);
    check_err(stat, __LINE__, __FILE__);
  }
  {
    // correlatlion functions
    size_t time_block = 4096 * 12 * 2 / ncells_;
    if (time_block < 1)
      time_block = 1;
    else if (time_block > get_n_frames())
      time_block = get_n_frames();
    size_t chunk[3] = {time_block, nrows_ / 2, ncols_};
    stat = nc_def_var_chunking(ncid_, corr_rho_id_, NC_CHUNKED, chunk);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_def_var_deflate(ncid_, corr_rho_id_, 0, 1, deflate_level);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_def_var_chunking(ncid_, corr_v_id_, NC_CHUNKED, chunk);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_def_var_deflate(ncid_, corr_v_id_, 0, 1, deflate_level);
    check_err(stat, __LINE__, __FILE__);
  }

  /* put global values */
  put_global_para(ncid_);

  /* leave define mode */
  stat = nc_enddef(ncid_);
  check_err(stat, __LINE__, __FILE__);
}

void SpatialCorr2dExporter::write_frame(int i_step, const float *c_rho,
                                        const float *c_v, const double *v_mean) {
  /* time step */
  {
  auto stat = nc_put_var1(ncid_, time_id_, time_idx_, &i_step);
  check_err(stat, __LINE__, __FILE__);
  }
  /* mean velocities */
  {
  size_t startset[2] = {time_idx_[0], 0};
  size_t countset[2] = {1, 2};
  auto stat = nc_put_vara(ncid_, v_mean_id_, startset, countset, &v_mean[0]);
  check_err(stat, __LINE__, __FILE__);
  }
  /* correlation functions */
  {
  size_t startset[3] = {time_idx_[0], 0, 0};
  size_t countset[3] = {1, nrows_ / 2, ncols_};
  auto stat = nc_put_vara(ncid_, corr_rho_id_, startset, countset, &c_rho[0]);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_put_vara(ncid_, corr_v_id_, startset, countset, &c_v[0]);
  check_err(stat, __LINE__, __FILE__);
  }
  ++time_idx_[0]; //! update time frame
}



