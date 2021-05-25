cmake_minimum_required(VERSION 3.16)
project(complete_matching)

set(CMAKE_CXX_STANDARD 20)
include_directories(include)

add_library(hashing STATIC include/hashing/hash_combine.h
                           source/hashing/hash_combine.cpp
                           source/hashing/hash_combine.cpp
                           include/hashing/pair_hash.h
                           include/hashing/vector_hash.h)
add_library(structures STATIC include/basics/trace.h
                              source/basics/trace.cpp
                              source/basics/AlterTrace.cpp
                              include/basics/AlterTrace.h)

add_executable(complete_matching main.cpp )
target_link_libraries(complete_matching hashing structures)