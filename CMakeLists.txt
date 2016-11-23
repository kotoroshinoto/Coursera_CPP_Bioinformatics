cmake_minimum_required(VERSION 3.6)
project(Coursera2)

set(CMAKE_CXX_STANDARD 11)
add_executable(translate translate.cpp)
add_executable(peptide_numseqs peptnumseq.cpp)
add_executable(frame_search frame_search.cpp)
add_executable(linear_spectrum linear_spectrum.cpp)
add_executable(cyclic_spectrum cyclic_spectrum.cpp)