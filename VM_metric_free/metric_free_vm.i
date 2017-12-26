%define DOCSTRING
""" Simulating metric-free vicsek model in 2D
"""
%enddef

%module (docstring=DOCSTRING) metric_free_vm
%{
    /* the resulting C file should be built as a python extension */
    #define SWIG_FILE_WITH_INIT
    /* Includes the header in the wrapper code */
    #include "metric_free_vm.h"

%}

%feature("autodoc", "1");

%include "metric_free_vm.h"