cmake_minimum_required(VERSION 2.8.12.2)
project(hpc C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wno-unknown-pragmas -Werror -pedantic-errors -std=c11")

add_executable(lab01 lab01.c)