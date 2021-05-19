cmake_minimum_required(VERSION 3.19)
project(lab5)

set(CMAKE_CXX_STANDARD 14)

find_package(MPI)
include_directories(SYSTEM ${MSMPI_INC})

add_executable(lab5 main.cpp)

target_link_libraries(${PROJECT_NAME} ${MPI_LIBRARIES})