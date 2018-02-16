cmake_minimum_required (VERSION 3.10.2)

MESSAGE (STATUS "Entering: " ${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory (mpi)
add_subdirectory (parafem_petsc)
add_subdirectory (choose_solvers)
