cmake_minimum_required(VERSION 3.3)
project(pagerank)
add_executable(pagerank pagerank.cxx)
SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3" )
