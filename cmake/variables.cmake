set ( BLAS_LIBRARIES
    "/usr/local/lib/libopenblas.a"
)
set ( BLAS_INCLUDE_DIR
    "/usr/local/include/"
)

set ( TPSM_LIBRARIES
    "/usr/local/lib/libtpsm.a"
)
set ( TPSM_INCLUDE_DIR
    "/usr/local/include/"
)

set ( LAPACKE_LIBRARIES
    "/usr/local/lib/liblapacke.a"
)
set ( LAPACKE_INCLUDE_DIR
    "/usr/local/include/"
)

set ( LAPACK_LIBRARIES
    "/usr/local/lib/liblapack.a"
)
set ( LAPACK_INCLUDE_DIR
    "/usr/local/include/"
)

include_directories(${BLAS_INCLUDE_DIR} ${TPSM_INCLUDE_DIR} ${LAPACKE_INCLUDE_DIR} ${LAPACK_INCLUDE_DIR})
