%define DOCSTRING
"Simulating 2D XY model"
%enddef

%module (docstring=DOCSTRING) XYmodel
%{
    /* the resulting C file should be built as a python extension */
    #define SWIG_FILE_WITH_INIT
    /* Includes the header in the wrapper code */
    #include "XYmodel.h"
%}

/* include the numpy typemaps */
// %ifdef _MSC_VER
%include "../VMpy/numpy.i"
// %else
// %include "numpy.i"
// %endif

/* need this for correct module initialization */
%init %{
  import_array();
%}

%feature("autodoc", "1");

%apply (double *IN_ARRAY1, int DIM1) {(double *_phi, int n_phi),
                                      (double *_psi, int n_psi)};

%apply (double *INPLACE_ARRAY1, int DIM1) {(double *phi_out, int size_phi)};

%rename (ini) my_ini;

%inline %{
  void my_ini(int _L, double _dt, double _eta, int seed,
              double *_phi, int n_phi, double _h, double *_psi, int n_psi) {
                ini(_L, _dt, _eta, seed, _phi, n_phi, _h, _psi);
              }
%}

%include "XYmodel.h"