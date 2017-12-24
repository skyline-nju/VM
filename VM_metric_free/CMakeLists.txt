# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( VM_metric_free )

cmake_minimum_required(VERSION 2.8.10)

find_package(CGAL QUIET)

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE "../../include")
  set(SRC_LIST vm_metric_free.cpp bird.cpp)
  add_executable(metric_free ${SRC_LIST})

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()