cmake_minimum_required(VERSION 3.12.3)
project(hpc C)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif ()

set(CMAKE_C_STANDARD 11)

set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3")

add_compile_options("-Wall")
add_compile_options("-Wno-unknown-pragmas")
add_compile_options("-Werror")
add_compile_options("-pedantic-errors")

add_executable(lab01 lab01.c)