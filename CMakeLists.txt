cmake_minimum_required(VERSION 2.8.12.2)
project(hpc C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_FLAGS "-Wall -Wno-unknown-pragmas -Werror -std=c11")

add_executable(hpc main.c)