cmake_minimum_required (VERSION 3.10.2)
MESSAGE (STATUS "Entering: " ${CMAKE_CURRENT_SOURCE_DIR})

file (GLOB ParaFEM_PETSc_SRC "*.F90")

add_library (parafem_petsc ${ParaFEM_PETSc_SRC}) 
add_dependencies(parafem_petsc mpi_module)
LinkFortranLibraries(parafem_petsc mpi_module)
