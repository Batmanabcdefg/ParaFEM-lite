#!/bin/bash

# This is just a quick and dirty little script that runs CMake, making some sensible-ish
# assumptions about your system setup

export FFLAGS="-O3 -march=native ${FFLAGS}"

cd build 
cmake ..
make -j8
