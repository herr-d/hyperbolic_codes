cmake_minimum_required (VERSION 3.9)
project(hyperbolic_codes)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE RelWithDebInfo)


set(CMAKE_CXX_FLAGS_DEBUG "-DDEBUG")

# ================ compiler
set(CMAKE_CXX_FLAGS "-std=c++11")


# ================ IMPORT THE AUTOTUNE LINK_LIBRARIES
# === IMPORTANT: AUTOTUNE MUST BE COMPLETELY COMPILED!!
# === CHANGE THE DIRECTORY VARIABLE TO THE AUTOTUNE FOLDER
set(AUTOTUNE_DIR "${PROJECT_SOURCE_DIR}/../autotune")
set(BLOSSOM_DIR "${PROJECT_SOURCE_DIR}/../autotune/blossomv")

# autotune include directory
include_directories(${AUTOTUNE_DIR}/libautotune/src/)

# add autotune link library
ADD_LIBRARY(autotune STATIC IMPORTED)
SET_TARGET_PROPERTIES(autotune PROPERTIES IMPORTED_LOCATION ${AUTOTUNE_DIR}/libautotune/lib/libautotune.a)

# add blossom object files
FILE(GLOB BLOSSOM_PM ${BLOSSOM_DIR}/PM*.o)
SET_SOURCE_FILES_PROPERTIES(
  ${BLOSSOM_PM}
  PROPERTIES
  EXTERNAL_OBJECT true
  GENERATED true
)
FILE(GLOB BLOSSOM_MISC ${BLOSSOM_DIR}/misc.o)
SET_SOURCE_FILES_PROPERTIES(
  ${BLOSSOM_MISC}
  PROPERTIES
  EXTERNAL_OBJECT true
  GENERATED true
)
FILE(GLOB BLOSSOM_MIN_COST ${BLOSSOM_DIR}/MinCost/MinCost.o)
SET_SOURCE_FILES_PROPERTIES(
  ${BLOSSOM_MIN_COST}
  PROPERTIES
  EXTERNAL_OBJECT true
  GENERATED true
)
FILE(GLOB BLOSSOM_GEOM ${BLOSSOM_DIR}/GEOM/*.o)
SET_SOURCE_FILES_PROPERTIES(
  ${BLOSSOM_GEOM}
  PROPERTIES
  EXTERNAL_OBJECT true
  GENERATED true
)

# define autotune error model folder (preprocessor variable)
add_definitions(-DEMS_DIR="${AUTOTUNE_DIR}/ex/ems/" )


# ================ include dirs
include_directories(include)

# ================ sub dirs
add_subdirectory(src)


# ================ testing
enable_testing()
add_subdirectory(Tests)
