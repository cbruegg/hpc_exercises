cmake_minimum_required(VERSION 3.12.3)
project(hpc C)

set(CMAKE_C_STANDARD 11)
add_compile_options("-Wall")
add_compile_options("-Wno-unknown-pragmas")
add_compile_options("-Werror")
add_compile_options("-pedantic-errors")

add_executable(lab01 lab01.c)