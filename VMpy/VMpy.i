%module VMpy
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

%apply (double* INPLACE_ARRAY1, int DIM1) {(double *x, int nBird),
                                           (double *y, int ny),
                                           (double *vx, int nvx),
                                           (double *vy, int nvy),
                                           (double *theta, int ncells)}

%rename (run) my_run;
%inline %{
    void my_run(int nstep,
                double *x, int nBird, double *y, int ny,
                double *vx, int nvx, double *vy, int nvy){
        run(nstep, x, y, vx, vy, nBird);
    }
%}

%rename (coarse_grain) my_coarse_grain;
%inline %{
    void my_coarse_grain(double l, double *theta, int ncells,
                         double *x, int nBird, double *y, int ny,
                         double *vx, int nvx, double *vy, int nvy) {
                             coarse_grain(l, theta, ncells, x, y, vx, vy, nBird);
                         }
%}

%include "VMpy.h"