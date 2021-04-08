cmake_minimum_required(VERSION 3.15)
project(lab1_1nd)

set(CMAKE_CXX_STANDARD 14)

find_package(MPI)
include_directories(SYSTEM ${MSMPI_INC})

add_executable(lab1_1nd main.cpp Matrix.h Constants.h)

target_link_libraries(${PROJECT_NAME} ${MPI_LIBRARIES})