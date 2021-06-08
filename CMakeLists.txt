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

add_library(graphs OBJECT include/graphs/TransitionGraph.h
                          source/graphs/adjacency_graph.cpp
                          include/graphs/adjacency_graph.h source/graphs/algorithms/nfa_to_dfa_weighted_labelled_automata.cpp include/graphs/algorithms/nfa_to_dfa_weighted_labelled_automata.h)

add_subdirectory(submodules/minauto)
include_directories(submodules/minauto/include/minauto)
include_directories(submodules/minauto/include)
include_directories(submodules/minauto/submodules/mpreal)
include_directories(submodules/minauto/submodules/ublas/include)
include_directories(submodules/minauto/submodules/boost_core/include)

add_executable(complete_matching main.cpp )
target_link_libraries(complete_matching equiv mpfr gmp hashing structures graphs)