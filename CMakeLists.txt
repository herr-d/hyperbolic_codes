cmake_minimum_required (VERSION 2.8)
project(hyperbolic_codes)

# ================ compiler
set(CMAKE_CXX_FLAGS "-std=c++11")

# ================ include dirs
include_directories(include)

# ================ sub dirs
add_subdirectory(graph_source)
add_subdirectory(src)
add_subdirectory(Tests)
