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

add_executable(Lab1
        src/lab1.cpp
        include/Timer.h
        include/Plotter.h
        include/IDataGenerator.h
        include/ISolver.h
        include/LabTask.h
        include/RandomSLAEGenerator.h
        include/ThomasSolver.h
        include/SLAE_ExactMethod_Task.h
        include/AsymptoticFuncs.h
        include/SLAE_IterMethod_Task.h
        include/JacobiSolver.h
        include/SeidelSolver.h)

target_link_libraries(Lab1 PRIVATE Labs)



