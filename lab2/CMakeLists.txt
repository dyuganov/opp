cmake_minimum_required(VERSION 3.17)
project(lab2)

set(CMAKE_CXX_STANDARD 14 )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -fopenmp")

add_compile_options(-mavx512f)

#add_executable(lab2 lab2.cpp)
#add_executable(lab2 lab2_vec.cpp)
add_executable(lab2 lab2_vec_parallel.cpp)
