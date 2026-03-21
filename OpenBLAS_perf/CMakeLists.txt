cmake_minimum_required(VERSION 3.10)
project(TRSM)

set(CMAKE_CXX_FLAGS "-O3 -march=native -ffast-math -fopenmp")

find_library(OPENBLAS_LIB openblas)

add_executable(benchmark
    benchmark/benchmark.cpp
    src/trsm_parallel.cpp
)

target_link_libraries(benchmark ${OPENBLAS_LIB} gomp)