
set(NUMA_FLAGS "-DNUMA_ALLOC")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -g -fPIC -O2 -fopenmp -D_GNU_SOURCE")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -lpthread -lnuma -DDLONG")

if (MERGE_DEFAULT)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DMERGE_DEFAULT ")
endif()

# when write graph features
if (WRITE_GRAPH)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DWRITE_GRAPH ")
endif()

include_directories("${CMAKE_SOURCE_DIR}/sparse/include")

include_directories("${CMAKE_SOURCE_DIR}/include")


set(LIB_BASE_SRCS
  src/base/SparseBase_config.c
  src/base/xerbla.c
  src/base/amd.c
  src/base/colamd.c
  )

set(LIB_CORE_SRCS
  src/core/SparseCore_common.c
  src/core/SparseCore_change_factor.c
  src/core/SparseCore_matrix_type.c
  src/core/SparseCore_matrixops.c
  src/core/SparseCore_check.c
  src/core/SparseCore_read_write.c
  src/core/SparseCore_norm.c
  # src/core/SparseCore_metis.c
  )

set(LIB_CHOL_SRCS
  src/chol/SparseChol_solve.c
  src/chol/SparseChol_super_solve.c
  src/chol/SparseChol_analyze.c
  src/chol/SparseChol_factorize.c
  src/chol/SparseChol_super_symbolic.c
  src/chol/SparseChol_super_numeric.c
  )

add_library(sparse STATIC
  ${LIB_BASE_SRCS} 
  ${LIB_CORE_SRCS} 
  ${LIB_CHOL_SRCS}
  ) 

target_link_libraries(sparse
  ${LAPACKE_LIBRARIES_DEP}
  ${CBLAS_LIBRARIES_DEP}
)
# Examples executables
if (WRITE_GRAPH)
set (EXAMPLES
  test/choltest.c
  test/cholwrite.c
  )
else()
set (EXAMPLES
  test/choltest.c
  )
endif()

foreach (_file ${EXAMPLES})
  get_filename_component(_name_we ${_file} NAME_WE)
  add_executable(sparse_${_name_we} ${_file})
  target_link_libraries(sparse_${_name_we} sparse  ${LAPACKE_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${TPSM_LIBRARIES} pthread gfortran numa m)
  # target_link_libraries(sparse_${_name_we} sparse  ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${TPSM_LIBRARIES} pthread gfortran numa m)
endforeach()

if (WRITE_GRAPH)
set (TEST_BASH 
    ./raw
    ./alldata.txt
    ./GCNdata_write.txt
    ./writeG.sh
    ./testChol.sh
)
else()
set (TEST_BASH 
    ./alldata.txt
    ./testChol.sh
)
endif()

foreach (_file ${TEST_BASH})
  file(COPY ${_file} DESTINATION ${CMAKE_BINARY_DIR}/test)
endforeach()