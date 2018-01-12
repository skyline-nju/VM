# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project( VM_metric_free )

cmake_minimum_required(VERSION 2.8.10)

find_package(CGAL QUIET)

option(LINK_STATIC_LIBS
       "Link the executable to static libraries"
       ON)
if (CGAL_FOUND)

  include( ${CGAL_USE_FILE} )

  #include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE "../../include")

  if (LINK_STATIC_LIBS)
    # link static libraries on unix.
    if (UNIX)
      set(my_flag "-static-libgcc -static-libstdc++ -Wl,-Bstatic -lmpfr -lgmp")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE} ${my_flag}")
      message(${CMAKE_CXX_FLAGS})
    endif(UNIX)
  endif(LINK_STATIC_LIBS)

  set(SRC_LIST main.cpp comn.cpp vm.cpp particle.cpp io_data.cpp)
  add_executable(metric_free ${SRC_LIST})

else()
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
endif()