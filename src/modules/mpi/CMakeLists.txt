cmake_minimum_required (VERSION 3.10.2)
MESSAGE (STATUS "Entering: " ${CMAKE_CURRENT_SOURCE_DIR})

file (GLOB MPI_MOD_SRC "*.f90")
file (GLOB SHARED_SRC "../shared/*.f90")

add_library (mpi_module ${MPI_MOD_SRC} ${SHARED_SRC}) 
