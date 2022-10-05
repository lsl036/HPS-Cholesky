cmake_minimum_required(VERSION 2.8)
project(hnuLapack)
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_COMPILER "gcc")
add_definitions("-O2")

set(USE_CONFIG_FILE ON)

message(STATUS "USE_CONFIG_FILE is ${USE_CONFIG_FILE}")
if(USE_CONFIG_FILE STREQUAL ON) 
    include(cmake/variables.cmake)
else()
    include(cmake/check.cmake)
endif()

set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/test)
add_subdirectory(${PROJECT_SOURCE_DIR}/sparse)

# add sparse libraries
# IF (ONLY_AMD STREQUAL ON)
#     ADD_DEFINITIONS(-DONLY_AMD)
# ENDIF()

# if (ONLY_METIS STREQUAL ON)
#     add_definitions(-DONLY_METIS)
# endif()

# if (ONLY_COLAMD STREQUAL ON)
#     add_definitions(-DONLY_COLAMD)
# endif()


