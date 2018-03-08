cmake_minimum_required (VERSION 3.7.0)
MESSAGE (STATUS "Entering: " ${CMAKE_CURRENT_SOURCE_DIR})

set (shared_modules shared_data.F90
                    plasticity.F90
                    setup.F90
                    linalg.F90
                    core.F90)

set (EXECUTABLES microFE)

foreach (EXECUTABLE ${EXECUTABLES})
  add_executable(${EXECUTABLE} ${EXECUTABLE}.F90 ${shared_modules})
  LinkFortranLibraries(${EXECUTABLE} parafem_petsc choose_solvers)
  LinkMPI(${EXECUTABLE})
  LinkPETSc(${EXECUTABLE})
  LinkLAPACK(${EXECUTABLE})
endforeach ()
