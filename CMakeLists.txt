cmake_minimum_required(VERSION 3.28)
project(numerical-methods-in-physics)

set(CMAKE_CXX_STANDARD 20)
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

function(add_tasks_from_glob glob_pattern)
    file(GLOB TASK_SOURCES "${glob_pattern}")

    foreach(task_src ${TASK_SOURCES})
        get_filename_component(task_name ${task_src} NAME_WE)
        add_executable(${task_name} ${task_src} ${LAB_COMMON_HEADERS})
        target_include_directories(${task_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
        target_link_libraries(${task_name} PRIVATE Labs)
    endforeach()
endfunction()

function(add_lab lab_number)
    set(is_special OFF)
    if(ARGC GREATER 1)
        set(is_special ${ARGV1})
    endif()

    if(is_special)
        set(pattern "${CMAKE_CURRENT_SOURCE_DIR}/src/Labs/Special/Lab${lab_number}/Task*.cpp")
    else()
        set(pattern "${CMAKE_CURRENT_SOURCE_DIR}/src/Labs/Lab${lab_number}/Task*.cpp")
    endif()

    add_tasks_from_glob("${pattern}")
endfunction()

function(special_add)
    add_tasks_from_glob(
            "${CMAKE_CURRENT_SOURCE_DIR}/src/Labs/Special/Additional/Task*.cpp"
    )
endfunction()

add_lab(6 ON)