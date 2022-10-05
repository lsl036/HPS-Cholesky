if(NOT DEFINED TPSM_PATH)
    find_package(PkgConfig)
    pkg_search_module(TPSM REQUIRED tpsm)
    if (TPSM_FOUND)
        message(STATUS "thread-pool : tpsm:${TPSM_VERSION}")
        include_directories( ${TPSM_INCLUDE_DIRS} )
        link_directories( ${TPSM_LIBRARY_DIRS} )
    else()
        message(FATAL_ERROR "can't not find the package tpsm, you can define TPSM_PATH point it")
    endif ()
else()
    include_directories(${TPSM_PATH}/include)
    link_directories( ${TPSM_PATH}/lib )
endif()
set(TPSM_LIBRARIES "tpsm")

if(NOT DEFINED BLAS_TYPE)
    set(BLAS_TYPE "openblas")
endif()

if(NOT DEFINED BLAS_PATH)
    find_package(PkgConfig)
    pkg_search_module(BLAS ${BLAS_TYPE})
    if (BLAS_FOUND)
        message(STATUS "blas : ${BLAS_TYPE}:${BLAS_VERSION}")
        include_directories( ${BLAS_INCLUDE_DIRS} )
        link_directories( ${BLAS_LIBRARY_DIRS} )
    else()
        find_path(BLAS_INCLUDE NAMES cblas.h PATHS /usr/include /usr/local/include /usr/local/${BLAS_TYPE}/include)    
        find_library(BLAS_LIB NAMES ${BLAS_TYPE} /usr/lib /usr/local/lib /usr/local/${BLAS_TYPE}/lib)
        message(STATUS "blas : ${BLAS_INCLUDE}/cblas.h and ${BLAS_LIB}")
        if (${BLAS_INCLUDE} STREQUAL "BLAS_INCLUDE-NOTFOUND" OR ${BLAS_LIB} STREQUAL "BLAS_LIB-NOTFOUND") 
            message(FATAL_ERROR "can't not find the package blas, you can define BLAS_PATH point it or define BLAS_TYPE to point which blas libraries to use")
        endif()
        include_directories(${BLAS_INCLUDE}) 
        set(BLAS_LIBRARIES ${BLAS_LIB})
    endif()
else()
    include_directories(${BLAS_PATH}/include) 
    link_directories( ${BLAS_PATH}/lib )
    set(BLAS_LIBRARIES ${BLAS_TYPE})
endif()    

if (NOT DEFINED LAPACKE_PATH)
    find_package(PkgConfig)
    pkg_search_module(LAPACKE lapacke)
    if (LAPACKE_FOUND) 
        message(STATUS "lapacke : lapacke:${LAPACKE_VERSION}")
        include_directories( ${LAPACKE_INCLUDE_DIRS} )
        link_directories( ${LAPACKE_LIBRARY_DIRS} )
    else()
        find_path(LAPACKE_INCLUDE NAMES lapacke.h PATHS /usr/include /usr/local/include /usr/local/lapack/include)    
        find_library(LAPACKE_LIB NAMES lapacke /usr/lib /usr/local/lib /usr/local/lapack/lib)
        find_library(LAPACK_LIB NAMES lapack /usr/lib /usr/local/lib /usr/local/lapack/lib)
        message(STATUS "lapacke : ${LAPACKE_INCLUDE}/lapacke.h and ${LAPACKE_LIB}")
        if (${LAPACKE_INCLUDE} STREQUAL "LAPACKE_INCLUDE-NOTFOUND" OR ${LAPACKE_LIB} STREQUAL "LAPACKE_LIB-NOTFOUND" OR ${LAPACK_LIB} STREQUAL "LAPACK_LIB-NOTFOUND") 
            message(FATAL_ERROR "can't not find the package lapacke, you can define LAPACKE_PATH to point it")
        endif()
        include_directories(${LAPACKE_INCLUDE})
        set(LAPACKE_LIBRARIES ${LAPACKE_LIB})
        set(LAPACK_LIBRARIES ${LAPACK_LIB})
    endif()
else()
    include_directories(${LAPACKE_PATH}/include)
    link_directories(${LAPACKE_PATH}/lib)
    set(LAPACKE_LIBRARIES "lapacke")
    set(LAPACK_LIBRARIES "lapack")
endif()


# if (NOT DEFINED METIS_PATH)
#     find_package(PkgConfig)
#     pkg_search_module(METIS metis)
#     if (METIS_FOUND) 
#         message(STATUS "metis : metis:${METIS_VERSION}")
#         include_directories( ${METIS_INCLUDE_DIRS} )
#         link_directories( ${METIS_LIBRARY_DIRS} )
#     else()
#         find_path(METIS_INCLUDE NAMES metis.h PATHS /usr/include /usr/local/include /usr/local/metis/include)    
#         find_library(METIS_LIB NAMES metis /usr/lib /usr/local/lib /usr/local/metis/lib)
#         message(STATUS "metis : ${METIS_INCLUDE}/metis.h and ${METIS_LIB}")  
#         if (${METIS_INCLUDE} STREQUAL "METIS_INCLUDE-NOTFOUND" OR ${METIS_LIB} STREQUAL "METIS_LIB-NOTFOUND") 
#             message(FATAL_ERROR "can't not find the package metis, you can define METIS_PATH to point it")
#         endif()
#         include_directories(${METIS_INCLUDE}) 
#         set(METIS_LIBRARIES ${METIS_LIB})   
#     endif()
# else()
#     include_directories(${METIS_PATH}/include)
#     link_directories(${METIS_PATH}/lib)
#     set(METIS_LIBRARIES "metis")
# endif()
