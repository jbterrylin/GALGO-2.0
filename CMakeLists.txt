# cmake -DCMAKE_BUILD_TYPE=Debug ..
# cmake -DCMAKE_BUILD_TYPE=Release ..
cmake_minimum_required(VERSION 3.0)
project(galgo)

# Variables
set(SOURCE_CORE ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(SOURCE_TEST ${CMAKE_CURRENT_SOURCE_DIR}/test)

# Compiler warning and optimization flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

# TODO - galgo Library - all headers only (*.hpp) for now
# file(GLOB GALGO_SOURCES
#   ${SOURCE_CORE}/*.hpp
# )
# add_library(galgo STATIC ${GALGO_SOURCES})

# Includes
include_directories(
  ${SOURCE_CORE}
)

# Test core
add_executable(testga ${SOURCE_CORE}/example.cpp)

# Test binairo
add_executable(testbinairo ${SOURCE_TEST}/Binairo/main.cpp)

