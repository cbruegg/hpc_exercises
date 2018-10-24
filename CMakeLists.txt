cmake_minimum_required(VERSION 3.10)
project(hpc C)

set(CMAKE_C_STANDARD 11)

add_compile_options("-Wall")
add_compile_options("-Wno-unknown-pragmas")
add_compile_options("-Werror")

add_executable(hpc main.c)