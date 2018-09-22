cmake_minimum_required(VERSION 3.12)
project(KD_Tree C)

set(CMAKE_C_STANDARD 99)

add_executable(KD_Tree main.c lib/utilities.c lib/utilities.h lib/kdTree.c lib/kdTree.h)