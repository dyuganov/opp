cmake_minimum_required(VERSION 3.15)
project(lab1)

set(CMAKE_CXX_STANDARD 14)

find_package(MPI)
include_directories(SYSTEM ${MSMPI_INC})

#add_executable(lab1 nonMPI.cpp)
#add_executable(lab1 MPI_v1.cpp Matrix.h)
#add_executable(lab1 MPI_v2.cpp Matrix.h)
add_executable(lab1 MPI_v2_2.cpp Matrix.h)

target_link_libraries(${PROJECT_NAME} ${MPI_LIBRARIES})
