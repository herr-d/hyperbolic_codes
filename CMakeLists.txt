cmake_minimum_required (VERSION 3.10)
project(hyperbolic_codes)


set(CMAKE_CXX_FLAGS_DEBUG "-DDEBUG")

# ================ compiler
set(CMAKE_CXX_FLAGS "-std=c++11")


# ================ IMPORT THE AUTOTUNE LINK_LIBRARIES
# === IMPORTANT: AUTOTUNE MUST BE COMPLETELY COMPILED!!
# === CHANGE THE DIRECTORY VARIABLE TO THE AUTOTUNE FOLDER
set(AUTOTUNE_DIR "${PROJECT_SOURCE_DIR}/../autotune")

# autotune include directory
include_directories(${AUTOTUNE_DIR}/libautotune/src/)

# add link library
ADD_LIBRARY(autotune STATIC IMPORTED)
SET_TARGET_PROPERTIES(autotune PROPERTIES IMPORTED_LOCATION ${AUTOTUNE_DIR}/libautotune/lib/libautotune.a)

# define autotune error model folder (preprocessor variable)
add_definitions(-DEMS_DIR="${AUTOTUNE_DIR}/ex/ems/" )


# ================ include dirs
include_directories(include)

# ================ sub dirs
add_subdirectory(src)


# ================ testing
enable_testing()
add_subdirectory(Tests)
