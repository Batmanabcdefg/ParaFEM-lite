cmake_minimum_required (VERSION 3.10.2)
MESSAGE (STATUS "Entering: " ${CMAKE_CURRENT_SOURCE_DIR})

file (GLOB choose_solvers_SRC "*.F90")

add_library (choose_solvers ${choose_solvers_SRC}) 
add_dependencies(choose_solvers mpi_module)
LinkFortranLibraries(choose_solvers mpi_module)
