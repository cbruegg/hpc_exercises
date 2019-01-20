cmake_minimum_required(VERSION 3.12.3)
project(hpc)

set(CMAKE_VERBOSE_MAKEFILE on)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif ()

set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 17)

add_compile_options("-Wall")
add_compile_options("-Wno-unknown-pragmas")
add_compile_options("-Werror")
add_compile_options("-pedantic-errors")
add_compile_options("-pthread")
add_compile_options("-fopenmp")

link_libraries("-fopenmp")

set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g ${COMPILE_OPTIONS}")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3 ${COMPILE_OPTIONS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g ${COMPILE_OPTIONS}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 ${COMPILE_OPTIONS}")

add_executable(lab01 lab01.c)
add_executable(lab02 lab02.cpp)
add_executable(lab03_1 lab03_1.cpp)
add_executable(lab03_2 lab03_2.cpp)
add_executable(lab04 lab04.cpp)
add_executable(lab05 lab05.cpp)
add_executable(lab06 lab06.cpp)
add_executable(lab07 lab07.cpp)
add_executable(lab08 lab08_src/RushHour.cpp lab08_src/Car.cpp lab08_src/State.cpp lab08_src/StateManager.cpp)
add_executable(lab09 lab09.cpp)
