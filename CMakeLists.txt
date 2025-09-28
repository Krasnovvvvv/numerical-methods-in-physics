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

file(GLOB_RECURSE LAB_COMMON_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp")

function(add_lab_task name)
    add_executable(${name} src/${name}.cpp ${LAB_COMMON_HEADERS})
    target_link_libraries(${name} PRIVATE Labs)
endfunction()

add_lab_task(Task1)



