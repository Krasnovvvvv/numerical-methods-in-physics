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

function(add_lab lab_number)
    set(is_special OFF)
    if(ARGC GREATER 1)
        set(is_special ${ARGV1})
    endif()

    if(is_special)
        file(GLOB TASK_SOURCES
                "${CMAKE_CURRENT_SOURCE_DIR}/src/Labs/Special/Lab${lab_number}/Task*.cpp"
        )
    else()
        file(GLOB TASK_SOURCES
                "${CMAKE_CURRENT_SOURCE_DIR}/src/Labs/Lab${lab_number}/Task*.cpp"
        )
    endif()

    foreach(task_src ${TASK_SOURCES})
        get_filename_component(task_name ${task_src} NAME_WE)
        add_executable(${task_name} ${task_src} ${LAB_COMMON_HEADERS})
        target_include_directories(${task_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
        target_link_libraries(${task_name} PRIVATE Labs)
    endforeach()
endfunction()


add_lab(6)


