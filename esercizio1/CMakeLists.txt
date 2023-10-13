cmake_minimum_required(VERSION 3.22)
project(SommaParallela)

find_package(MPI REQUIRED)

add_executable(sommaParallela sommaPar.c)
target_link_libraries(sommaParallela m)