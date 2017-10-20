%define DOCSTRING
"Simulating Vicsek model in 2D."
%enddef

%module (docstring=DOCSTRING) VMpy
%{
    /* the resulting C file should be built as a python extension */
    #define SWIG_FILE_WITH_INIT
    /* Includes the header in the wrapper code */
    #include "VMpy.h"
%}

/* include the numpy typemaps */
%include "numpy.i"

/* need this for correct module initialization */
%init %{
    import_array();
%}

%feature("autodoc", "1");

%apply (double* INPLACE_ARRAY1, int DIM1) {(double *x, int nBird),
                                           (double *y, int ny),
                                           (double *vx, int nvx),
                                           (double *vy, int nvy),
                                           (double *theta, int ncells),
                                           (double *svx, int nsvx),
                                           (double *svy, int nsvy)}

%apply (int* INPLACE_ARRAY1, int DIM1) {(int *num, int len_num)}

%rename (ini) my_ini;
%inline %{
  void my_ini(double *x, int nBird, double *y, int ny,
              double *vx, int nvx, double *vy, int nvy,
              int seed, double v0, double _eta, double Lx, double Ly) {
                ini(x, y, vx, vy, nBird, seed, v0, _eta, Lx, Ly);
  }
%}

%rename (get_snap) my_get_snap;
%inline %{
  void my_get_snap(double *x, int nBird, double *y, int ny,
                   double *vx, int nvx, double *vy, int nvy) {
                     get_snap(x, y, vx, vy, nBird);
  }
%}

%rename (get_coarse_grained_snap) my_get_coarse_grained_snap;
%inline %{
  void my_get_coarse_grained_snap(int *num, int len_num, double *svx, int nsvx,
                                  double *svy, int nsvy, double l) {
                                    get_coarse_grained_snap(num, svx, svy, len_num, l);
  }
%}

%include "VMpy.h"