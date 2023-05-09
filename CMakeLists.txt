cmake_minimum_required(VERSION 2.8)
project(hnuLapack)
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_COMPILER "gcc")
add_definitions("-O2")

set(MERGE_DEFAULT ON)
set(USE_CONFIG_FILE OFF)
set(WRITE_GRAPH OFF)    # collect the graph's properties and generate raw dataset

message(STATUS "MERGE_DEFAULT is ${MERGE_DEFAULT}")
message(STATUS "USE_CONFIG_FILE is ${USE_CONFIG_FILE}")
message(STATUS "WRITE_GRAPH is ${WRITE_GRAPH}")

if(USE_CONFIG_FILE STREQUAL ON) 
    include(cmake/variables.cmake)
else()
    include(cmake/check.cmake)
endif()

set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/test)
add_subdirectory(${PROJECT_SOURCE_DIR}/sparse)



