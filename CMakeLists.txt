cmake_minimum_required(VERSION 3.28)
project(numerical-methods-in-physics)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 CONFIG REQUIRED)
find_package(Matplot++ CONFIG REQUIRED)

add_library(Labs INTERFACE)

target_link_libraries(Labs INTERFACE
        Eigen3::Eigen
        Matplot++::cimg
        Matplot++::matplot)

add_executable(Lab1 lab1.cpp
        lab1.cpp
        Timer.h
        Plotter.h)
target_link_libraries(Lab1 PRIVATE Labs)



