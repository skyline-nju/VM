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
%include "../VMpy/numpy.i"

/* need this for correct module initialization */
%init %{
  import_array();
%}

%feature("autodoc", "1");

%apply (double *IN_ARRAY1, int DIM1) {(double *_phi, int n_phi)};

%apply (double *INPLACE_ARRAY1, int DIM1) {(double *phi_out, int size_phi)};

%include "XYmodel.h"