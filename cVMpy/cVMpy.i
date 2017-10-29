%define DOCSTRING
"""Simulating 2D XY model with continuum time."""
%enddef

%module (docstring=DOCSTRING) cVMpy
%{
    /* the resulting C file should be built as a python extension */
    #define SWIG_FILE_WITH_INIT
    /* Includes the header in the wrapper code */
    #include "cVMpy.h"
%}

/* include the numpy typemaps */
%include "../VMpy/numpy.i"

/* need this for correct module initialization */
%init %{
  import_array();
%}

%feature("autodoc", "1");

%apply (double *IN_ARRAY1, int DIM1) {(double *x_in, int nx_in),
                                     (double *y_in, int ny_in),
                                     (double *theta_in, int ntheta_in)}

%apply (double *INPLACE_ARRAY1, int DIM1) {(double *x_out, int nx_out),
                                          (double *y_out, int ny_out),
                                          (double *theta_out, int ntheta_out)}

%rename (ini) my_ini;
%rename (get_snap) my_get_snap;

%inline %{
  void my_ini(double *x_in, int nx_in, double *y_in, int ny_in,
              double *theta_in, int ntheta_in, double _eta, double Lx, double Ly,
              int seed, double dt) {
                ini(x_in, y_in, theta_in, nx_in, _eta, Lx, Ly, seed, dt);
              }
  void my_get_snap(double *x_out, int nx_out, double *y_out, int ny_out,
                   double *theta_out, int ntheta_out){
      get_snap(x_out, y_out, theta_out, nx_out);
    }
%}

%include "cVMpy.h"
