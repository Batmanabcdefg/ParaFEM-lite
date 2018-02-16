cmake_minimum_required (VERSION 3.10.2)
project (ParaFEM)
enable_language(Fortran)

set (CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

list (APPEND CMAKE_MODULE_PATH ${ParaFEM_SOURCE_DIR}/cmake-modules)

add_definitions(-DREVISION=\"PT001\")

find_package (PETSc REQUIRED)
find_package (MPI REQUIRED)
find_package(LAPACK REQUIRED)

include_directories (${PETSC_INCLUDES} ${MPI_Fortran_INCLUDE_PATH})

function(LinkFortranLibraries Target)
  target_link_libraries(${Target} ${ARGN})
  foreach(Lib IN LISTS ARGN)
    target_include_directories(${Target} PUBLIC $<TARGET_PROPERTY:${Lib},Fortran_MODULE_DIRECTORY>)
    MESSAGE (STATUS "Linking: " $<TARGET_PROPERTY:${Lib},Fortran_MODULE_DIRECTORY>)
  endforeach()
endfunction()

function(LinkPETSc Target)
  target_link_libraries(${Target} ${MPI_Fortran_LIBRARIES})
endfunction()

function(LinkMPI Target)
  target_link_libraries(${Target} ${PETSC_LIBRARIES})
endfunction()

function(LinkLapack Target)
  target_link_libraries(${Target} ${LAPACK_LIBRARIES})
endfunction()

MESSAGE (STATUS ${PETSC_LIBRARIES})

add_subdirectory (${ParaFEM_SOURCE_DIR}/src)
