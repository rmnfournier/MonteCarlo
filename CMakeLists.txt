cmake_minimum_required(VERSION 3.10)
project(MonteCarlo)
set(CMAKE_C_COMPILER /usr/bin/mpicc)
set(CMAKE_CXX_COMPILER /usr/bin/mpic++)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-g -pg -O3 -ftree-vectorize -funroll-loops")
add_executable(MonteCarlo main.cpp MonteCarlo.cpp MonteCarlo.hh Interacting_Harmonic_Oscillator.cpp Interacting_Harmonic_Oscillator.hh)