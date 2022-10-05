#ifndef _models_h_
#define _models_h_

#include "plu_common.h"
#include "plu_solver.h"
#include "plu_flops.h"

#ifndef _kernels_trace_h_
#define _kernels_trace_h_

/**
 * @brief Main stop enum event for all the events in traces
 */
#define PluKernelStop 0

/**
 * @brief List of the Level 0 events that may be traced in PLU
 *
 * This is only the high level steps.
 */
typedef enum plu_ktype0_e {
    PluKernelLvl0Facto,
    PluKernelLvl0Solve,
    PluKernelLvl0Diag,
    PluKernelLvl0Nbr
} plu_ktype0_t;

/**
 * @brief List of the Level 1 events that may be traced in PLU
 *
 * This is the main information that traces all the major kernels during the
 * factorization step.
 */
typedef enum plu_ktype_e {
    PluKernelGETRF,        /**< LU diagonal block kernel             */
    PluKernelSCALOCblk,    /**< Scaling out-of-place of a panel      */
    PluKernelSCALOBlok,    /**< Scaling out-of-place of a block      */
    PluKernelTRSMCblk1d,   /**< TRSM applied to a panel in 1d layout */
    PluKernelTRSMCblk2d,   /**< TRSM applied to a panel in 2d layout */
    PluKernelTRSMCblkLR,   /**< TRSM applied to a panel in low-rank  */
    PluKernelTRSMBlok2d,   /**< TRSM applied to a block in 2d layout */
    PluKernelTRSMBlokLR,   /**< TRSM applied to a block in low-rank  */
    PluKernelGEMMCblk1d1d, /**< GEMM applied from a panel in 1d layout to a panel in 1d layout */
    PluKernelGEMMCblk1d2d, /**< GEMM applied from a panel in 1d layout to a panel in 2d layout */
    PluKernelGEMMCblk2d2d, /**< GEMM applied from a panel in 2d layout to a panel in 2d layout */
    PluKernelGEMMCblkFRLR, /**< GEMM applied from a panel in full-rank to a panel in low-rank  */
    PluKernelGEMMCblkLRLR, /**< GEMM applied from a panel in low-rank to a panel in low-rank   */
    PluKernelGEMMBlok2d2d, /**< GEMM applied from a block in 2d layout to a block in 2d layout */
    PluKernelGEMMBlokLRLR, /**< GEMM applied from a block in low-rank to a block in low-rank   */
    PluKernelGEADDCblkFRFR, /**< GEADD applied from a panel in full-rank to a panel in full-rank  */
    PluKernelGEADDCblkFRLR, /**< GEADD applied from a panel in full-rank to a panel in low-rank  */
    PluKernelGEADDCblkLRLR, /**< GEADD applied from a panel in low-rank to a panel in low-rank   */
    PluKernelLvl1Nbr
} plu_ktype_t;

/**
 * @brief List of the Level 2 events that may be traced in PLU
 *
 * This is the low-level information that traces all the individual calls to
 * blas/lapack routines in the code. It is used to compute the number of flops
 * in low-rank compression, and to distinguish the amount of flops spent in each
 * part of the low-rank updates.
 *
 */
typedef enum plu_ktype2_e {

    /* General kernels: similar in low-rank and dense */
    PluKernelLvl2GETRF,             /**< LU diagonal block kernel           */
    PluKernelLvl2HETRF,             /**< LDLh diagonal block kernel         */
    PluKernelLvl2POTRF,             /**< Cholesky diagonal block kernel     */
    PluKernelLvl2PXTRF,             /**< Complex LL^t diagonal block kernel */
    PluKernelLvl2SYTRF,             /**< LDLt diagonal block kernel         */

    /* Solve operations */
    PluKernelLvl2_FR_TRSM,
    PluKernelLvl2_LR_TRSM,

    /* Update operations */
    PluKernelLvl2_FR_GEMM,

    /* Formation (and application) of A * B */
    PluKernelLvl2_LR_FRFR2FR,
    PluKernelLvl2_LR_FRLR2FR,
    PluKernelLvl2_LR_LRFR2FR,
    PluKernelLvl2_LR_LRLR2FR,
    PluKernelLvl2_LR_FRFR2LR,
    PluKernelLvl2_LR_FRLR2LR,
    PluKernelLvl2_LR_LRFR2LR,
    PluKernelLvl2_LR_LRLR2LR,
    PluKernelLvl2_LR_FRFR2null,
    PluKernelLvl2_LR_FRLR2null,
    PluKernelLvl2_LR_LRFR2null,
    PluKernelLvl2_LR_LRLR2null,

    /* Compression kernels */
    PluKernelLvl2_LR_init_compress,
    PluKernelLvl2_LR_add2C_uncompress,
    PluKernelLvl2_LR_add2C_recompress,
    PluKernelLvl2_LR_add2C_updateCfr,
    PluKernelLvl2_LR_add2C_orthou,
    PluKernelLvl2_LR_add2C_rradd_orthogonalize, /**<< CGS, partialQR or fullQR */
    PluKernelLvl2_LR_add2C_rradd_recompression,
    PluKernelLvl2_LR_add2C_rradd_computeNewU,

    PluKernelLvl2Nbr
} plu_ktype2_t;

/**
 * @brief Total number of kernel events
 */
#define PluKernelsNbr (PluKernelLvl0Nbr + PluKernelLvl1Nbr + PluKernelLvl2Nbr)

/**
 * @brief Global array to store the number of flops executed per kernel
 */
extern volatile double kernels_flops[PluKernelLvl1Nbr];

/**
 * @brief Lock to accumulate flops
 */
extern plu_atomic_lock_t lock_flops;

/**
 * @brief Overall number of flops
 */
extern double overall_flops[3];

#if defined(PLU_WITH_EZTRACE)

#include "eztrace_module/kernels_ev_codes.h"

/**
 * @brief Define the level traced by the EZTrace module
 */
extern int plu_eztrace_level;

#else

static inline void kernel_trace_start_lvl0     ( plu_ktype0_t ktype )  { (void)ktype; }
static inline void kernel_trace_stop_lvl0      ( double flops )           { (void)flops; }
static inline void kernel_trace_start_lvl2     ( plu_ktype2_t ktype )  { (void)ktype; }
static inline void kernel_trace_stop_lvl2      ( double flops )           { (void)flops; }
static inline void kernel_trace_stop_lvl2_rank ( double flops, int rank ) { (void)flops; (void)rank; }

#endif

#if defined(PLU_GENERATE_MODEL)

/**
 * @brief Structure to store information linked to a kernel in order to generate
 * the cost model
 */
typedef struct plu_model_entry_s {
    plu_ktype_t ktype; /**< The type of the kernel             */
    int m;                /**< The first diemension of the kernel */
    int n;                /**< The second dimension of the kernel if present, 0 otherwise */
    int k;                /**< The third dimension of the kernel, 0 otherwise             */
    double time;          /**< The time spent in the kernel (s)                           */
} plu_model_entry_t;

extern plu_model_entry_t *model_entries;     /**< Array to all entries                 */
extern volatile int32_t      model_entries_nbr; /**< Index of the last entry in the array */
extern int32_t               model_size;        /**< Size of the model_entries array      */

#endif

void   kernelsTraceStart( const plu_data_t *plu_data );
double kernelsTraceStop(  const plu_data_t *plu_data );

/**
 *******************************************************************************
 *
 * @brief Start the trace of a single kernel
 *
 *******************************************************************************
 *
 * @param[in] ktype
 *          Type of the kernel starting that need to be traced.
 *          With EZTrace mode, this call is empty if the environment variable
 *          PLU_EZTRACE_LEVEL is different from 1.
 *
 *******************************************************************************
 *
 * @return the starting time if PLU_GENERATE_MODEL is enabled, 0. otherwise.
 *
 *******************************************************************************/
static inline double
kernel_trace_start( plu_ktype_t ktype )
{
    double time = 0.;

#if defined(PLU_WITH_EZTRACE)

    if (plu_eztrace_level == 1) {
        EZTRACE_EVENT_PACKED_0( KERNELS_LVL1_CODE(ktype) );
    }

#endif

#if defined(PLU_GENERATE_MODEL)

    time = clockGet();

#endif

    (void)ktype;
    return time;
}

/**
 *******************************************************************************
 *
 * @brief Stop the trace of a single kernel
 *
 *******************************************************************************
 *
 * @param[in] ktype
 *          Type of the kernel starting that need to be traced.
 *          With EZTrace mode, this call is empty if the environment variable
 *          PLU_EZTRACE_LEVEL is different from 1.
 *
 * @param[in] m
 *          The m parameter of the kernel (used by xxTRF, TRSM, and GEMM)
 *
 * @param[in] n
 *          The n parameter of the kernel (used by TRSM, and GEMM)
 *
 * @param[in] k
 *          The k parameter of the kernel (used by GEMM)
 *
 * @param[in] flops
 *          The number of flops of the kernel
 *
 * @param[in] starttime
 *          The stating time of the kernel. Used only if PLU_GENERATE_MODEL
 *          is enabled.
 *
 *******************************************************************************/
static inline void
kernel_trace_stop( int8_t inlast, plu_ktype_t ktype, int m, int n, int k, double flops, double starttime )
{

#if defined(PLU_WITH_EZTRACE)

    if (plu_eztrace_level == 1) {
        EZTRACE_EVENT_PACKED_1( KERNELS_CODE( PluKernelStop ), flops );
    }

#endif

#if defined(PLU_GENERATE_MODEL)

    {
        double  time = clockGet() - starttime;
        int32_t index = plu_atomic_inc_32b( &model_entries_nbr );

        if ( index < model_size ) {
            model_entries[index].ktype = ktype;
            model_entries[index].m     = m;
            model_entries[index].n     = n;
            model_entries[index].k     = k;
            model_entries[index].time  = time;
        }
        else {
            fprintf(stderr, "WARNING: too many kernels to log %d\n", index);
        }
    }

#endif

    plu_atomic_lock( &lock_flops );
    overall_flops[inlast] += flops;
    plu_atomic_unlock( &lock_flops );

    (void)ktype;
    (void)m;
    (void)n;
    (void)k;
    (void)flops;
    (void)starttime;
    return;
}

#endif /* _kernels_trace_h_ */


/**
 * @brief Model structure to store the coefficients table and its name
 */
typedef struct plu_model_s {
    char  *name;                                     /**< Name of the computational unit considered by the model */
    double coefficients[4][PluKernelLvl1Nbr][8];  /**< Coefficients table of the model                        */
} plu_model_t;

/**
 * @brief Return the time in s of factorization kernels using a single size parameter.
 * @param[in] coefs The coefficients array to use in the formula
 * @param[in] N The size parameter
 * @return The estimated time in s.
 */
static inline double
modelsGetCost1Param( const double *coefs, Sparse_long N )
{
    /* a3 * N^3 + a2 * N^2 + a1 * N + a0 */
    double time = ((coefs[3] * N + coefs[2]) * N + coefs[1]) * N + coefs[0];
    return (time < 0.) ? 0. : time;
}

/**
 * @brief Return the time in s of TRSM kernels using two size parameters.
 * @param[in] coefs The coefficients array to use in the formula
 * @param[in] M The first size parameter (number of rows)
 * @param[in] N The second size parameter (number of columns, size of the triangular matrix)
 * @return The estimated time in s.
 */
static inline double
modelsGetCost2Param( const double *coefs, Sparse_long M, Sparse_long N )
{
    /* a5 * M*N^2 + a4 * M*N + a3 * N^2 + a2 * M + a1 * N + a0 */
    double time = ((coefs[5] * (double)M + coefs[3]) * (double)N + coefs[4] * (double)M + coefs[1]) * (double)N + coefs[2] * (double)M + coefs[0];
    return (time < 0.) ? 0. : time;
}

/**
 * @brief Return the time in s of GEMM kernels using three size parameters.
 * @param[in] coefs The coefficients array to use in the formula
 * @param[in] M The number of rows in the C matrix
 * @param[in] N The number of columns in the C matrix
 * @param[in] K The third dimension of the matrix product
 * @return The estimated time in s.
 */
static inline double
modelsGetCost3Param( const double *coefs, Sparse_long M, Sparse_long N, Sparse_long K )
{
    /* a7 * M * N * K + a6 * M * K + a5 * K * N + a4 * M * N + a3 * M + a2 * N + a1 * K + a0 */
    double time = (coefs[7] * (double)M * (double)N * (double)K +
                   coefs[6] * (double)M * (double)K +
                   coefs[5] * (double)K * (double)N +
                   coefs[4] * (double)M * (double)N +
                   coefs[3] * (double)M +
                   coefs[2] * (double)N +
                   coefs[1] * (double)K +
                   coefs[0]);
    return (time < 0.) ? 0. : time;
}

plu_model_t *pluModelsNew();
void pluModelsFree( plu_model_t *model );
void pluModelsLoad( plu_data_t *plu_data );

#endif /* _models_h_ */

#ifndef _sopalin_data_h_
#define _sopalin_data_h_

struct sopalin_data_s {
    SolverMatrix *solvmtx;
    double      (*cpu_coefs)[PluKernelLvl1Nbr][8];
    double      (*gpu_coefs)[PluKernelLvl1Nbr][8];
};
typedef struct sopalin_data_s sopalin_data_t;

void sopalin_ztrsm( plu_data_t *plu_data, int side, int uplo, int trans, int diag, sopalin_data_t *sopalin_data, int nrhs, plu_complex64_t *b, int ldb );
void sopalin_ctrsm( plu_data_t *plu_data, int side, int uplo, int trans, int diag, sopalin_data_t *sopalin_data, int nrhs, plu_complex32_t *b, int ldb );
void sopalin_dtrsm( plu_data_t *plu_data, int side, int uplo, int trans, int diag, sopalin_data_t *sopalin_data, int nrhs, double *b, int ldb );
void sopalin_strsm( plu_data_t *plu_data, int side, int uplo, int trans, int diag, sopalin_data_t *sopalin_data, int nrhs, float *b, int ldb );

void sopalin_zdiag( plu_data_t *plu_data, sopalin_data_t *sopalin_data, int nrhs, plu_complex64_t *b, int ldb );
void sopalin_cdiag( plu_data_t *plu_data, sopalin_data_t *sopalin_data, int nrhs, plu_complex32_t *b, int ldb );
void sopalin_ddiag( plu_data_t *plu_data, sopalin_data_t *sopalin_data, int nrhs, double *b, int ldb );
void sopalin_sdiag( plu_data_t *plu_data, sopalin_data_t *sopalin_data, int nrhs, float *b, int ldb );

void sopalin_zgetrf( plu_data_t *plu_data, sopalin_data_t *sopalin_data );
void sopalin_cgetrf( plu_data_t *plu_data, sopalin_data_t *sopalin_data );
void sopalin_dgetrf( plu_data_t *plu_data, sopalin_data_t *sopalin_data );
void sopalin_sgetrf( plu_data_t *plu_data, sopalin_data_t *sopalin_data );

void sopalin_zhetrf( plu_data_t *plu_data, sopalin_data_t *sopalin_data );
void sopalin_chetrf( plu_data_t *plu_data, sopalin_data_t *sopalin_data );

void sopalin_zpotrf( plu_data_t *plu_data, sopalin_data_t *sopalin_data );
void sopalin_cpotrf( plu_data_t *plu_data, sopalin_data_t *sopalin_data );
void sopalin_dpotrf( plu_data_t *plu_data, sopalin_data_t *sopalin_data );
void sopalin_spotrf( plu_data_t *plu_data, sopalin_data_t *sopalin_data );

void sopalin_zpxtrf( plu_data_t *plu_data, sopalin_data_t *sopalin_data );
void sopalin_cpxtrf( plu_data_t *plu_data, sopalin_data_t *sopalin_data );

void sopalin_zsytrf( plu_data_t *plu_data, sopalin_data_t *sopalin_data );
void sopalin_csytrf( plu_data_t *plu_data, sopalin_data_t *sopalin_data );
void sopalin_dsytrf( plu_data_t *plu_data, sopalin_data_t *sopalin_data );
void sopalin_ssytrf( plu_data_t *plu_data, sopalin_data_t *sopalin_data );

#endif /* _sopalin_data_h_ */

#ifndef _coeftab_d_h_
#define _coeftab_d_h_

Sparse_long coeftab_dcompress  ( SolverMatrix *solvmtx );
void         coeftab_duncompress( SolverMatrix *solvmtx );
void         coeftab_dmemory    ( SolverMatrix *solvmtx );

void coeftab_dgetschur( const SolverMatrix *solvmtx,
                        double *S, Sparse_long lds );

void coeftab_dgetdiag( const SolverMatrix *solvmtx,
                       double *D, Sparse_long incD );

void coeftab_ddump    ( plu_data_t      *plu_data,
                        const SolverMatrix *solvmtx,
                        const char         *filename );
int  coeftab_ddiff    ( plu_coefside_t   side,
                        const SolverMatrix *solvA,
                        SolverMatrix       *solvB );

#endif /* _coeftab_d_h_ */

#ifndef _coeftab_h_
#define _coeftab_h_

void coeftabInit( plu_data_t     *plu_data,
                  plu_coefside_t  side );
void coeftabExit( SolverMatrix      *solvmtx );

Sparse_long coeftabCompress( plu_data_t *plu_data );

/**
 * @brief Type of the memory gain functions
 */
typedef void (*coeftab_fct_memory_t)( SolverMatrix * );

/**
 * @brief List of functions to compute the memory gain in low-rank per precision.
 */
coeftab_fct_memory_t coeftabMemory[4];

#endif /* _coeftab_h_ */


#ifndef _d_nan_check_h_
#define _d_nan_check_h_

#ifndef DOXYGEN_SHOULD_SKIP_THIS
// #endif /* DOXYGEN_SHOULD_SKIP_THIS */

#if defined(PLU_DEBUG_LR_NANCHECK)
#define LAPACKE_dlacpy_work LAPACKE_dlacpy
#define LAPACKE_dlaset_work LAPACKE_dlaset

#define LAPACKE_dormlq_work( _layout_, _side_, _trans_, _m_, _n_, _k_, _a_, _lda_, _tau_, _c_, _ldc_, _w_, _ldw_ ) \
    LAPACKE_dormlq( _layout_, _side_, _trans_, _m_, _n_, _k_, _a_, _lda_, _tau_, _c_, _ldc_ )
#define LAPACKE_dormqr_work( _layout_, _side_, _trans_, _m_, _n_, _k_, _a_, _lda_, _tau_, _c_, _ldc_, _w_, _ldw_ ) \
    LAPACKE_dormqr( _layout_, _side_, _trans_, _m_, _n_, _k_, _a_, _lda_, _tau_, _c_, _ldc_ )

#define LAPACKE_dgeqrf_work( _layout_, _m_, _n_, _a_, _lda_, _tau_, _w_, _ldw_ ) \
    LAPACKE_dgeqrf( _layout_, _m_, _n_, _a_, _lda_, _tau_ )
#define LAPACKE_dgelqf_work( _layout_, _m_, _n_, _a_, _lda_, _tau_, _w_, _ldw_ ) \
    LAPACKE_dgelqf( _layout_, _m_, _n_, _a_, _lda_, _tau_ )

#if defined(PRECISION_z) || defined(PRECISION_c)
#define MYLAPACKE_dgesvd_work( _layout_, _jobu_, jobv_, _m_, _n_, _a_, _lda_, _s_, _u_, _ldu_, _v_, _ldv_, _w_, _ldw_, _rw_ ) \
    LAPACKE_dgesvd( _layout_, _jobu_, jobv_, _m_, _n_, _a_, _lda_, _s_, _u_, _ldu_, _v_, _ldv_, (double*)(_w_) )
#else
#define MYLAPACKE_dgesvd_work( _layout_, _jobu_, jobv_, _m_, _n_, _a_, _lda_, _s_, _u_, _ldu_, _v_, _ldv_, _w_, _ldw_, _rw_ ) \
    LAPACKE_dgesvd( _layout_, _jobu_, jobv_, _m_, _n_, _a_, _lda_, _s_, _u_, _ldu_, _v_, _ldv_, (double*)(_w_) )
#endif

#else

#if defined(PRECISION_z) || defined(PRECISION_c)
#define MYLAPACKE_dgesvd_work( _layout_, _jobu_, jobv_, _m_, _n_, _a_, _lda_, _s_, _u_, _ldu_, _v_, _ldv_, _w_, _ldw_, _rw_ ) \
    LAPACKE_dgesvd_work( _layout_, _jobu_, jobv_, _m_, _n_, _a_, _lda_, _s_, _u_, _ldu_, _v_, _ldv_, _w_, _ldw_, _rw_ )
#else
#define MYLAPACKE_dgesvd_work( _layout_, _jobu_, jobv_, _m_, _n_, _a_, _lda_, _s_, _u_, _ldu_, _v_, _ldv_, _w_, _ldw_, _rw_ ) \
    LAPACKE_dgesvd_work( _layout_, _jobu_, jobv_, _m_, _n_, _a_, _lda_, _s_, _u_, _ldu_, _v_, _ldv_, _w_, _ldw_ )
#endif

#endif /* defined(PLU_DEBUG_LR_NANCHECK) */

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#endif /* _d_nan_check_h_ */



#ifndef _d_refine_functions_h_
#define _d_refine_functions_h_

struct d_solver
{
    Sparse_long    (* getN   )   (plu_data_t *);
    plu_fixdbl_t (* getEps )   (plu_data_t *);
    Sparse_long    (* getImax)   (plu_data_t *);
    Sparse_long    (* getRestart)(plu_data_t *);

    void* (*malloc)(size_t);
    void  (*free)(void*);

    void   (*output_oneiter)(double, double, double, Sparse_long);
    void   (*output_final)( plu_data_t *, double, Sparse_long,
                            double, void*, double*);

    void   (*scal)( plu_data_t *, Sparse_long, double, double * );
    double (*dot) ( plu_data_t *, Sparse_long, const double *, const double * );
    void   (*copy)( plu_data_t *, Sparse_long, const double *, double * );
    void   (*axpy)( plu_data_t *, Sparse_long, double, const double *, double *);
    void   (*spmv)( const plu_data_t *, plu_trans_t, double, const double *, double, double * );
    void   (*spsv)( plu_data_t *, double * );
    double (*norm)( plu_data_t *, Sparse_long, const double * );
    void   (*gemv)( plu_data_t *, Sparse_long, Sparse_long,
                    double, const double *, Sparse_long,
                    const double *, double, double *);
};

void d_refine_init(struct d_solver *, plu_data_t*);

Sparse_long d_gmres_smp   ( plu_data_t *plu_data, void *x, void *b );
Sparse_long d_grad_smp    ( plu_data_t *plu_data, void *x, void *b );
Sparse_long d_pivot_smp   ( plu_data_t *plu_data, void *x, void *b );
Sparse_long d_bicgstab_smp( plu_data_t *plu_data, void *x, void *b );

#endif /* _d_refine_functions_h_ */


#ifndef _plu_dcores_h_
#define _plu_dcores_h_

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define plu_cblk_lock( cblk_ )    plu_atomic_lock( &((cblk_)->lock) )
#define plu_cblk_unlock( cblk_ )  plu_atomic_unlock( &((cblk_)->lock) )
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 * @addtogroup kernel_blas_lapack
 * @{
 *    This module contains all the BLAS and LAPACK-like kernels that are working
 *    on lapack layout matrices.
 *
 *    @name PluDouble BLAS kernels
 *    @{
 */
void core_dplrnt( int m, int n, double *A, int lda,
                  int gM, int m0, int n0, unsigned long long int seed );
void core_dgetmo( int m, int n,
                  const double *A, int lda,
                  double *B, int ldb );
int  core_dgeadd( plu_trans_t trans, Sparse_long M, Sparse_long N,
                  double alpha, const double *A, Sparse_long LDA,
                  double beta,        double *B, Sparse_long LDB );
int  core_dgemdm( plu_trans_t transA, plu_trans_t transB, int M, int N, int K,
                  double  alpha, const double *A, int LDA,
                  const double *B, int LDB,
                  double  beta, double *C, int LDC,
                  const double *D, int incD,
                  double *WORK, int LWORK );
int  core_dpqrcp( double tol, Sparse_long maxrank, int full_update, Sparse_long nb,
                  Sparse_long m, Sparse_long n,
                  double *A, Sparse_long lda,
                  Sparse_long *jpvt, double *tau,
                  double *work, Sparse_long lwork,  double *rwork );
int  core_drqrcp( double tol, Sparse_long maxrank, int refine, Sparse_long nb,
                  Sparse_long m, Sparse_long n,
                  double *A, Sparse_long lda,
                  Sparse_long *jpvt, double *tau,
                  double *work, Sparse_long lwork,  double *rwork );
int  core_drqrrt( double tol, Sparse_long maxrank, Sparse_long nb,
                  Sparse_long m, Sparse_long n,
                  double *A, Sparse_long lda, double *tau,
                  double *B, Sparse_long ldb, double *tau_b,
                  double *work, Sparse_long lwork,  double normA );
int  core_dtqrcp( double tol, Sparse_long maxrank, int unused, Sparse_long nb,
                  Sparse_long m, Sparse_long n,
                  double *A, Sparse_long lda,
                  Sparse_long *jpvt, double *tau,
                  double *work, Sparse_long lwork,  double *rwork );
int  core_dtradd( plu_uplo_t uplo, plu_trans_t trans, Sparse_long M, Sparse_long N,
                  double alpha, const double *A, Sparse_long LDA,
                  double beta,        double *B, Sparse_long LDB);
int  core_dscalo( plu_trans_t trans, Sparse_long M, Sparse_long N,
                  const double *A, Sparse_long lda,
                  const double *D, Sparse_long ldd,
                  double *B, Sparse_long ldb );

/**
 *    @}
 *    @name PluDouble Othogonalization kernels for low-rank updates
 *    @{
 */
plu_fixdbl_t
core_dlrorthu_fullqr( Sparse_long M,  Sparse_long N, Sparse_long rank,
                      double *U, Sparse_long ldu,
                      double *V, Sparse_long ldv );
plu_fixdbl_t
core_dlrorthu_partialqr( Sparse_long M,  Sparse_long N,
                         Sparse_long r1, Sparse_long *r2ptr,
                         Sparse_long offx, Sparse_long offy,
                         double *U, Sparse_long ldu,
                         double *V, Sparse_long ldv );
plu_fixdbl_t
core_dlrorthu_cgs( Sparse_long M1,  Sparse_long N1,
                   Sparse_long M2,  Sparse_long N2,
                   Sparse_long r1, Sparse_long *r2ptr,
                   Sparse_long offx, Sparse_long offy,
                   double *U, Sparse_long ldu,
                   double *V, Sparse_long ldv );

/**
 *    @}
 *    @name PluDouble LAPACK kernels
 *    @{
 */
void core_dpotrfsp( Sparse_long n, double *A, Sparse_long lda,
                    Sparse_long *nbpivot, double criterion );
void core_dpotrfsp( Sparse_long n, double *A, Sparse_long lda,
                    Sparse_long *nbpivot, double criterion );
void core_dgetrfsp( Sparse_long n, double *A, Sparse_long lda,
                    Sparse_long *nbpivot, double criterion );
void core_dsytrfsp( Sparse_long n, double *A, Sparse_long lda,
                    Sparse_long *nbpivot, double criterion );
void core_dsytrfsp( Sparse_long n, double *A, Sparse_long lda,
                    Sparse_long *nbpivot, double criterion );

/**
 *     @}
 * @}
 *
 * @addtogroup kernel_fact
 * @{
 *    This module contains all the kernel working at the solver matrix structure
 *    level for the numerical factorization step.
 *
 *    @name PluDouble cblk-BLAS CPU kernels
 *    @{
 */

int  cpucblk_dgeaddsp1d( const SolverCblk *cblk1, SolverCblk *cblk2,
                         const double *L1, double *L2,
                         const double *U1, double *U2 );

void cpucblk_dgemmsp( plu_coefside_t sideA, plu_coefside_t sideB, plu_trans_t trans,
                      const SolverCblk *cblk, const SolverBlok *blok, SolverCblk *fcblk,
                      const double *A, const double *B, double *C,
                      double *work, Sparse_long lwork, const plu_lr_t *lowrank );
void cpucblk_dtrsmsp( plu_coefside_t coef, plu_side_t side, plu_uplo_t uplo,
                      plu_trans_t trans, plu_diag_t diag, SolverCblk *cblk,
                      const double *A, double *C,
                      SolverMatrix *solvmtx );
void cpucblk_dscalo ( plu_trans_t trans, SolverCblk *cblk, double *LD );

void cpublok_dgemmsp( plu_coefside_t sideA, plu_coefside_t sideB, plu_trans_t trans,
                      const SolverCblk *cblk, SolverCblk *fcblk,
                      Sparse_long blok_mk, Sparse_long blok_nk, Sparse_long blok_mn,
                      const double *A, const double *B, double *C,
                      const plu_lr_t *lowrank );
void cpublok_dtrsmsp( plu_coefside_t coef, plu_side_t side, plu_uplo_t uplo,
                      plu_trans_t trans, plu_diag_t diag,
                      SolverCblk *cblk, Sparse_long blok_m,
                      const double *A, double *C,
                      const plu_lr_t *lowrank );
void cpublok_dscalo ( plu_trans_t trans,
                      SolverCblk *cblk, Sparse_long blok_m,
                      const double *A, const double *D, double *B );

/**
 *    @}
 *    @name PluDouble cblk LU kernels
 *    @{
 */
int cpucblk_dgetrfsp1d_getrf( SolverMatrix *solvmtx, SolverCblk *cblk,
                              double *L, double *U );
int cpucblk_dgetrfsp1d_panel( SolverMatrix *solvmtx, SolverCblk *cblk,
                              double *L, double *U );
int cpucblk_dgetrfsp1d      ( SolverMatrix *solvmtx, SolverCblk *cblk,
                              double *work, Sparse_long lwork );

/**
 *    @}
 *    @name PluDouble initialization and additionnal routines
 *    @{
 */
void cpucblk_dalloc   ( plu_coefside_t    side,
                        SolverCblk          *cblk );
void cpucblk_dfillin  ( plu_coefside_t    side,
                        const SolverMatrix  *solvmtx,
                        const plu_bcsc_t *bcsc,
                        Sparse_long         itercblk );
void cpucblk_dinit    ( plu_coefside_t    side,
                        const SolverMatrix  *solvmtx,
                        const plu_bcsc_t *bcsc,
                        Sparse_long         itercblk,
                        const char          *directory );
void cpucblk_dgetschur( const SolverCblk    *cblk,
                        int                  upper_part,
                        double  *S,
                        Sparse_long         lds );
void cpucblk_ddump    ( plu_coefside_t    side,
                        const SolverCblk    *cblk,
                        FILE                *stream );
int  cpucblk_ddiff    ( plu_coefside_t    side,
                        const SolverCblk    *cblkA,
                        SolverCblk          *cblkB );
void cpucblk_dadd     ( plu_coefside_t    side,
                        double               alpha,
                        const SolverCblk    *cblkA,
                        SolverCblk          *cblkB,
                        const plu_lr_t   *lowrank );
/**
 *    @name PluDouble compression/uncompression routines
 */
Sparse_long cpucblk_dcompress( const SolverMatrix *solvmtx,
                                plu_coefside_t   side,
                                SolverCblk         *cblk );
void         cpucblk_duncompress( plu_coefside_t side,
                                  SolverCblk       *cblk );
void         cpucblk_dmemory    ( plu_coefside_t  side,
                                  SolverMatrix      *solvmtx,
                                  SolverCblk        *cblk,
                                  Sparse_long      *gain );

/**
 *    This module contains all the kernel working on the solver matrix structure
 *    for the solve step.
 */

void solve_blok_dtrsm( plu_coefside_t coefside, plu_side_t side, plu_uplo_t uplo,
                       plu_trans_t trans, plu_diag_t diag, const SolverCblk *cblk,
                       int nrhs, double *b, int ldb );
void solve_blok_dgemm( plu_coefside_t coefside, plu_side_t side, plu_trans_t trans,
                       Sparse_long nrhs, const SolverCblk *cblk, const SolverBlok *blok,
                       SolverCblk *fcbk, const double *B, Sparse_long ldb,
                       double *C, Sparse_long ldc );

void solve_cblk_dtrsmsp_forward( plu_solv_mode_t mode, plu_side_t side, plu_uplo_t uplo,
                                 plu_trans_t trans, plu_diag_t diag,
                                 const SolverMatrix *datacode, const SolverCblk *cblk,
                                 int nrhs, double *b, int ldb );
void solve_cblk_dtrsmsp_backward( plu_solv_mode_t mode, plu_side_t side, plu_uplo_t uplo,
                                  plu_trans_t trans, plu_diag_t diag,
                                  const SolverMatrix *datacode, const SolverCblk *cblk,
                                  int nrhs, double *b, int ldb );

void solve_cblk_ddiag( const SolverCblk   *cblk,
                       int                 nrhs,
                       double *b,
                       int                 ldb,
                       double *work );

/**
 *    This module contains the three terms update functions for the LDL^t and
 *    LDL^h factorizations.
 */
void core_dsytrfsp1d_gemm( const SolverCblk *cblk, const SolverBlok *blok, SolverCblk *fcblk,
                           const double *L, double *C,
                           double *work );
void core_dsytrfsp1d_gemm( const SolverCblk *cblk, const SolverBlok *blok, SolverCblk *fcblk,
                           const double *L, double *C,
                           double *work );

#endif /* _plu_dcores_h_ */