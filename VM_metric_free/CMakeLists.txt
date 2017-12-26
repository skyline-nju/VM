# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project( VM_metric_free )

cmake_minimum_required(VERSION 2.8.10)

find_package(CGAL QUIET)

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE "../../include")

  # link static libraries on unix.
  if (UNIX)
    set(my_flag "-static-libgcc -static-libstdc++ -Wl,-Bstatic -lmpfr -lgmp")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE} ${my_flag}")
    message(${CMAKE_CXX_FLAGS})
  endif(UNIX)

  set(SRC_LIST main.cpp bird.cpp metric_free_vm.cpp comn.cpp)
  add_executable(metric_free ${SRC_LIST})

  # CAGL-SWIG bindings
  # find_package(SWIG REQUIRED)
  # include( ${SWIG_USE_FILE} )

  # find_package(PythonLibs)
  # include_directories( ${PYTHON_INCLUDE_PATH})

  # include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

  # set(CMAKE_SWIG_FLAGS "")

  # set_source_files_properties(metric_free_vm.i PROPERTIES CPLUSPLUS ON)
  # set_source_files_properties(metric_free_vm.i PROPERTIES SWIG_FLAGS "-includeall")
  # SWIG_ADD_LIBRARY(metric_free_vm LANGUAGE python  SOURCES metric_free_vm.i metric_free_vm.cpp bird.cpp)
  # SWIG_LINK_LIBRARIES(metric_free_vm ${PYTHON_LIBRARIES})

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()