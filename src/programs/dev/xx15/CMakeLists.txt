cmake_minimum_required (VERSION 3.7.0)
MESSAGE (STATUS "Entering: " ${CMAKE_CURRENT_SOURCE_DIR})

file (GLOB ParaFEM_PETSc_SRC "*.F90")

set (EXECUTABLES_xx15 xx15_elastic 
                      xx15_quadric_linear_hardening
                      xx15_vm_nonlinear_hardening)

foreach (EXECUTABLE ${EXECUTABLES_xx15})
  add_executable(${EXECUTABLE} ${EXECUTABLE}.F90 plasticity_xx15.F90)
  LinkFortranLibraries(${EXECUTABLE} parafem_petsc choose_solvers)
  LinkMPI(${EXECUTABLE})
  LinkPETSc(${EXECUTABLE})
  LinkLAPACK(${EXECUTABLE})
endforeach ()
