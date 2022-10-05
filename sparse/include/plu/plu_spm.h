
#ifndef _spm_h_
#define _spm_h_

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include "SparseCore.h"

#ifndef _spm_config_h_
#define _spm_config_h_

#define SPM_VERSION_MAJOR 
#define SPM_VERSION_MINOR 
#define SPM_VERSION_MICRO 

/* system */
#define HAVE_CLOCK_GETTIME
#define HAVE_ASPRINTF
#define HAVE_VASPRINTF
#define HAVE_STDARG_H
#define HAVE_UNISTD_H
#define HAVE_VA_COPY
/* #undef HAVE_UNDERSCORE_VA_COPY */
#define HAVE_GETOPT_LONG
#define HAVE_GETRUSAGE
#define HAVE_GETOPT_H
#define HAVE_ERRNO_H
#define HAVE_STDDEF_H
#define HAVE_LIMITS_H
#define HAVE_STRING_H
#define HAVE_COMPLEX_H
/* #undef HAVE_FALLTHROUGH */

/* Datatypes used */
#define SPM_INT64

/* LAPACKE */
/* #undef LAPACKE_WITH_LASCL */

#if defined(HAVE_FALLTHROUGH)
#define spm_attr_fallthrough __attribute__((fallthrough))
#else
#define spm_attr_fallthrough do {} while(0)
#endif

#if defined(WIN32) || defined(_WIN32)
#define SPM_OS_WINDOWS 1
#endif

/*
 * BEGIN_C_DECLS should be used at the beginning of your declarations,
 * so that C++ compilers don't mangle their names.  Use END_C_DECLS at
 * the end of C declarations.
 */
#undef BEGIN_C_DECLS
#undef END_C_DECLS
#if defined(c_plusplus) || defined(__cplusplus)
# define BEGIN_C_DECLS extern "C" {
# define END_C_DECLS }
#else
#define BEGIN_C_DECLS          /* empty */
#define END_C_DECLS            /* empty */
#endif

#endif /* _spm_config_h_ */

#ifndef _spm_api_h_
#define _spm_api_h_

/********************************************************************
 * CBLAS value address
 */
#ifndef CBLAS_SADDR
#define CBLAS_SADDR( a_ ) (&(a_))
#endif

/**
 * @brief Verbose modes
 */
typedef enum spm_verbose_e {
    SpmVerboseNot = 0, /**< Nothing  */
    SpmVerboseNo  = 1, /**< Default  */
    SpmVerboseYes = 2  /**< Extended */
} spm_verbose_t;

/**
 * @brief Error codes
 */
typedef enum spm_error_e {
    SPM_SUCCESS            = 0,  /**< No error                     */
    SPM_ERR_UNKNOWN        = 1,  /**< Unknown error                */
    SPM_ERR_ALLOC          = 2,  /**< Allocation error             */
    SPM_ERR_NOTIMPLEMENTED = 3,  /**< Not implemented feature      */
    SPM_ERR_OUTOFMEMORY    = 4,  /**< Not enough memory            */
    SPM_ERR_THREAD         = 5,  /**< Error with threads           */
    SPM_ERR_INTERNAL       = 6,  /**< Internal error               */
    SPM_ERR_BADPARAMETER   = 7,  /**< Bad parameters given         */
    SPM_ERR_FILE           = 8,  /**< Error in In/Out operations   */
    SPM_ERR_INTEGER_TYPE   = 9,  /**< Error with integer types     */
    SPM_ERR_IO             = 10, /**< Error with input/output      */
    SPM_ERR_MPI            = 11  /**< Error with MPI calls         */
} spm_error_t;

/**
 * @brief How to generate RHS
 */
typedef enum spm_rhstype_e {
    SpmRhsOne,
    SpmRhsI,
    SpmRhsRndX,
    SpmRhsRndB
} spm_rhstype_t;

/**
 * @brief Upper/Lower part
 */
typedef enum spm_uplo_e {
    SpmUpper      = 121, /**< Use lower triangle of A */
    SpmLower      = 122, /**< Use upper triangle of A */
    SpmUpperLower = 123  /**< Use the full A          */
} spm_uplo_t;

/**
 * @brief Diagonal
 */
typedef enum spm_diag_e {
    SpmNonUnit = 131, /**< Diagonal is non unitary */
    SpmUnit    = 132  /**< Diagonal is unitary     */
} spm_diag_t;

/**
 * @brief Side of the operation
 */
typedef enum spm_side_e {
    SpmLeft  = 141, /**< Apply operator on the left  */
    SpmRight = 142  /**< Apply operator on the right */
} spm_side_t;

/**
 * @brief Norms
 */
typedef enum spm_normtype_e {
    SpmOneNorm       = 171, /**< One norm:       max_j( sum_i( |a_{ij}| ) )   */
    SpmFrobeniusNorm = 174, /**< Frobenius norm: sqrt( sum_{i,j} (a_{ij}^2) ) */
    SpmInfNorm       = 175, /**< Inifinite norm: max_i( sum_j( |a_{ij}| ) )   */
    SpmMaxNorm       = 177  /**< Inifinite norm: max_{i,j}( | a_{ij} | )      */
} spm_normtype_t;

/**
 * @brief Direction
 */
typedef enum spm_dir_e {
    SpmDirForward  = 391, /**< Forward direction   */
    SpmDirBackward = 392, /**< Backward direction  */
} spm_dir_t;

#endif /* _spm_api_h_ */

#ifndef _spm_datatypes_h_
#define _spm_datatypes_h_

#include <inttypes.h>

static inline Sparse_long spm_imin( Sparse_long a, Sparse_long b ) {
    return ( a < b ) ? a : b;
}

static inline Sparse_long spm_imax( Sparse_long a, Sparse_long b ) {
    return ( a > b ) ? a : b;
}

static inline Sparse_long spm_iceil( Sparse_long a, Sparse_long b ) {
    return ( a + b - 1 ) / b;
}

/** ****************************************************************************
 * Double that are not converted through precision generator functions
 **/
typedef double spm_fixdbl_t;

static inline size_t
spm_size_of(spm_coeftype_t type)
{
    switch(type) {
    case SpmDouble:    return   sizeof(double);
    default:
        fprintf(stderr, "spm_size_of: invalid type parameter\n");
        assert(0);
        return sizeof(double);
    }
}

#endif /* _spm_datatypes_h_ */

/**
 * @name SPM basic subroutines
 */
void spmAlloc( sparse_csc *spm );
void spmExit ( sparse_csc *spm );

sparse_csc *spmCopy( const sparse_csc *spm );
void        spmBase( sparse_csc *spm, int baseval );
Sparse_long   spmFindBase( const sparse_csc *spm );
int         spmConvert( int ofmttype, sparse_csc *ospm );
void        spmUpdateComputedFields( sparse_csc *spm );
void        spmGenFakeValues( sparse_csc *spm );

/**
 * @name SPM BLAS subroutines
 */
double spmNorm( spm_normtype_t ntype, const sparse_csc *spm );
int    spmMatVec( spm_trans_t trans, double alpha, const sparse_csc *spm, const void *x, double beta, void *y );
int    spmMatMat( spm_trans_t trans, Sparse_long n,
                  double alpha, const sparse_csc *A,
                                const void *B, Sparse_long ldb,
                  double beta,        void *C, Sparse_long ldc );
void   spmScalMatrix( double alpha, sparse_csc *spm );
void   spmScalVector( spm_coeftype_t flt, double alpha, Sparse_long n, void *x, Sparse_long incx );

/**
 * @name SPM subroutines to check format
 */
int       spmSort( sparse_csc *spm );
Sparse_long spmMergeDuplicate( sparse_csc *spm );
Sparse_long spmSymmetrize( sparse_csc *spm );
int       spmCheckAndCorrect( const sparse_csc *spm_in, sparse_csc *spm_out );

/**
 * @name SPM subroutines to check factorization/solve
 */
int spmGenRHS( spm_rhstype_t     type,
               Sparse_long         nrhs,
               const sparse_csc *spm,
               void * const      x,
               Sparse_long         ldx,
               void * const      b,
               Sparse_long         ldb );
int spmCheckAxb( double            eps,
                 Sparse_long         nrhs,
                 const sparse_csc *spm,
                 void * const      x0,
                 Sparse_long         ldx0,
                 void * const      b,
                 Sparse_long         ldb,
                 const void *      x,
                 Sparse_long         ldx );

/**
 * @name SPM subroutines to manipulate integers arrays
 */
Sparse_long *spmIntConvert( Sparse_long n, int *input );
void       spmIntSort1Asc1( void *const pbase, const Sparse_long n );
void       spmIntSort2Asc1( void *const pbase, const Sparse_long n );
void       spmIntSort2Asc2( void *const pbase, const Sparse_long n );
void       spmIntMSortIntAsc( void **const pbase, const Sparse_long n );

/**
 * @name SPM IO subroutines
 */
int spmLoad( sparse_csc *spm, FILE *infile );
int spmSave( const sparse_csc *spm, FILE *outfile );

/**
 * @name SPM driver
 */
int spmReadDriver( const char  *filename,
                   sparse_csc  *spm );

/**
 * @name SPM debug subroutines
 */
void *      spm2Dense( const sparse_csc *spm );
void        spmPrint( const sparse_csc *spm, FILE *f );
void        spmPrintInfo( const sparse_csc *spm, FILE *f );
sparse_csc *spmExpand( const sparse_csc *spm );
sparse_csc *spmDofExtend( const sparse_csc *spm, const int type, const int dof );

/**
 * @name SPM dev printing subroutines
 */

/**
 * @brief Subroutines to print one element of an spm structure
 * @param[in] f Pointer to the file
 * @param[in] i Row index of the element
 * @param[in] j Column index of the element
 * @param[in] A Value of the element A|i,j]
 * @details Double real case
 */
static inline void
d_spmPrintElt( FILE *f, Sparse_long i, Sparse_long j, double A )
{
    fprintf( f, "%ld %ld %e\n", (long)i, (long)j, A );
}

/**
 * @copydoc z_spmPrintElt
 * @details Pattern case
 *
 * @remark: uses a macro to avoid accessing A that would generate segfault.
 */
#define p_spmPrintElt( f, i, j, A )                                                                \
    {                                                                                              \
        fprintf( f, "%ld %ld\n", (long)( i ), (long)( j ) );                                       \
    }

/**
 * @}
 */
#endif /* _spm_h_ */

#ifndef _spm_common_h_
#define _spm_common_h_

#include <unistd.h>
#include <assert.h>
#include <errno.h>
#include <inttypes.h>
#include <limits.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>

#if defined(SPM_OS_WINDOWS)
#include <windows.h>
#endif

/********************************************************************
 * Errors functions
 */
#if defined(__GNUC__)
static inline void spm_print_error  ( const char *fmt, ...) __attribute__((format(printf,1,2)));
static inline void spm_print_warning( const char *fmt, ...) __attribute__((format(printf,1,2)));
#endif

static inline void
spm_print_error( const char *fmt, ... )
{
    va_list arglist;
    va_start(arglist, fmt);
    vfprintf(stderr, fmt, arglist);
    va_end(arglist);
}

static inline void
spm_print_warning( const char *fmt, ... )
{
    va_list arglist;
    va_start(arglist, fmt);
    fprintf(stderr, "WARNING: ");
    vfprintf(stderr, fmt, arglist);
    va_end(arglist);
}

/********************************************************************
 * CBLAS value address
 */
#ifndef CBLAS_SADDR
#define CBLAS_SADDR( a_ ) (&(a_))
#endif

/*
 * Get environment variable
 */
#if defined(SPM_OS_WINDOWS)

static inline int
spm_setenv( const char *var, const char *value, int overwrite ) {
    return !(SetEnvironmentVariable( var, value ));
}

static inline char *
spm_getenv( const char *var ) {
    char *str;
    int len = 512;
    int rc;
    str = (char*)malloc(len * sizeof(char));
    rc = GetEnvironmentVariable(var, str, len);
    if (rc == 0) {
        free(str);
        str = NULL;
    }
    return str;
}

static inline void
spm_cleanenv( char *str ) {
    if (str != NULL) free(str);
}

#else /* Other OS systems */

static inline int
spm_setenv( const char *var, const char *value, int overwrite ) {
    return setenv( var, value, overwrite );
}

static inline char *
spm_getenv( const char *var ) {
    return getenv( var );
}

static inline void
spm_cleanenv( char *str ) {
    (void)str;
}

#endif

#endif /* _spm_common_h_ */



#ifndef _spm_drivers_h_
#define _spm_drivers_h_

void convertArrayToDouble(    Sparse_long n, const double *A, void **B );

int readMM   ( const char *filename, sparse_csc *spm );
int readDMM  ( const char *filename, sparse_csc *spm );

#endif /* _spm_drivers_h_ */


#ifndef _p_spm_h_
#define _p_spm_h_

/**
 * Integer routines
 */
void p_spmIntFltSortAsc(void ** const pbase, const Sparse_long n);
void p_spmIntIntFltSortAsc(void ** const pbase, const Sparse_long n);

/**
 * Conversion routines
 */
int p_spmConvertCSC2CSR( sparse_csc *spm );
int p_spmConvertCSC2IJV( sparse_csc *spm );
int p_spmConvertCSR2CSC( sparse_csc *spm );
int p_spmConvertCSR2IJV( sparse_csc *spm );
int p_spmConvertIJV2CSC( sparse_csc *spm );
int p_spmConvertIJV2CSR( sparse_csc *spm );

int *p_spm2dense( const sparse_csc *spm );

/**
 * Matrix-Vector and matrix-matrix product routines
 */
int spm_pspmv( spm_trans_t            trans,
               int        alpha,
               const sparse_csc      *A,
               const int *x,
               Sparse_long              incx,
               int        beta,
               int       *y,
               Sparse_long              incy );
int spm_pspmm( spm_side_t             side,
               spm_trans_t            transA,
               spm_trans_t            transB,
               Sparse_long              K,
               int        alpha,
               const sparse_csc      *A,
               const int *B,
               Sparse_long              ldb,
               int        beta,
               int       *C,
               Sparse_long              ldc );

/**
 * Norm computation routines
 */
int p_spmNorm( int ntype, const sparse_csc *spm );

/**
 * Extra routines
 */
void      p_spmSort( sparse_csc *spm );
Sparse_long p_spmMergeDuplicate( sparse_csc *spm );
Sparse_long p_spmSymmetrize( sparse_csc *spm );

int p_spmGenRHS(spm_rhstype_t type, int nrhs, const sparse_csc *spm, void *x, int ldx, void *b, int ldb );
int p_spmCheckAxb( spm_fixdbl_t eps, int nrhs, const sparse_csc *spm, void *x0, int ldx0, void *b, int ldb, const void *x, int ldx );

/**
 * Output routines
 */
void p_spmDensePrint( FILE *f, Sparse_long m, Sparse_long n, const int *A, Sparse_long lda );
void p_spmPrint( FILE *f, const sparse_csc *spm );

sparse_csc *p_spmExpand(const sparse_csc *spm);
void        p_spmDofExtend(sparse_csc *spm);
void        p_spmScal( const int alpha, sparse_csc *spm );


#endif /* _p_spm_h_ */


#ifndef _d_spm_h_
#define _d_spm_h_

/**
 * Integer routines
 */
void d_spmIntFltSortAsc(void ** const pbase, const Sparse_long n);
void d_spmIntIntFltSortAsc(void ** const pbase, const Sparse_long n);

/**
 * Conversion routines
 */
int d_spmConvertCSC2CSR( sparse_csc *spm );
int d_spmConvertCSC2IJV( sparse_csc *spm );
int d_spmConvertCSR2CSC( sparse_csc *spm );
int d_spmConvertCSR2IJV( sparse_csc *spm );
int d_spmConvertIJV2CSC( sparse_csc *spm );
int d_spmConvertIJV2CSR( sparse_csc *spm );

double *d_spm2dense( const sparse_csc *spm );

/**
 * Matrix-Vector and matrix-matrix product routines
 */
int spm_dspmv( spm_trans_t            trans,
               double        alpha,
               const sparse_csc      *A,
               const double *x,
               Sparse_long              incx,
               double        beta,
               double       *y,
               Sparse_long              incy );
int spm_dspmm( spm_side_t             side,
               spm_trans_t            transA,
               spm_trans_t            transB,
               Sparse_long              K,
               double        alpha,
               const sparse_csc      *A,
               const double *B,
               Sparse_long              ldb,
               double        beta,
               double       *C,
               Sparse_long              ldc );

/**
 * Norm computation routines
 */
double d_spmNorm( int ntype, const sparse_csc *spm );

/**
 * Extra routines
 */
void      d_spmSort( sparse_csc *spm );
Sparse_long d_spmMergeDuplicate( sparse_csc *spm );
Sparse_long d_spmSymmetrize( sparse_csc *spm );

int d_spmGenRHS(spm_rhstype_t type, int nrhs, const sparse_csc *spm, void *x, int ldx, void *b, int ldb );
int d_spmCheckAxb( spm_fixdbl_t eps, int nrhs, const sparse_csc *spm, void *x0, int ldx0, void *b, int ldb, const void *x, int ldx );

/**
 * Output routines
 */
void d_spmDensePrint( FILE *f, Sparse_long m, Sparse_long n, const double *A, Sparse_long lda );
void d_spmPrint( FILE *f, const sparse_csc *spm );

sparse_csc *d_spmExpand(const sparse_csc *spm);
void        d_spmDofExtend(sparse_csc *spm);
void        d_spmScal( const double alpha, sparse_csc *spm );


#endif /* _d_spm_h_ */

/*
*   Matrix Market I/O library for ANSI C
*
*   See http://math.nist.gov/MatrixMarket for details.
*
*
*/

#ifndef _mmio_h_
#define _mmio_h_

#define MM_MAX_LINE_LENGTH 1025
#define MatrixMarketBanner "%%MatrixMarket"
#define MM_MAX_TOKEN_LENGTH 64

typedef char MM_typecode[4];

char *mm_typecode_to_str(MM_typecode matcode);

int mm_read_banner(FILE *f, MM_typecode *matcode);
int mm_read_mtx_crd_size(FILE *f, int *M, int *N, int *nz);
int mm_read_mtx_array_size(FILE *f, int *M, int *N);

int mm_write_banner(FILE *f, MM_typecode matcode);
int mm_write_mtx_crd_size(FILE *f, int M, int N, int nz);
int mm_write_mtx_array_size(FILE *f, int M, int N);


/********************* MM_typecode query fucntions ***************************/

#define mm_is_matrix(typecode)     ((typecode)[0]=='M')

#define mm_is_sparse(typecode)     ((typecode)[1]=='C')
#define mm_is_coordinate(typecode) ((typecode)[1]=='C')
#define mm_is_dense(typecode)      ((typecode)[1]=='A')
#define mm_is_array(typecode)      ((typecode)[1]=='A')

#define mm_is_real(typecode)       ((typecode)[2]=='R')
#define mm_is_pattern(typecode)    ((typecode)[2]=='P')
#define mm_is_integer(typecode)    ((typecode)[2]=='I')

#define mm_is_symmetric(typecode)  ((typecode)[3]=='S')
#define mm_is_general(typecode)    ((typecode)[3]=='G')

int mm_is_valid(MM_typecode matcode);		/* too complex for a macro */


/********************* MM_typecode modify fucntions ***************************/

#define mm_set_matrix(typecode)      ((*typecode)[0]='M')
#define mm_set_coordinate(typecode)  ((*typecode)[1]='C')
#define mm_set_array(typecode)       ((*typecode)[1]='A')
#define mm_set_dense(typecode)       mm_set_array(typecode)
#define mm_set_sparse(typecode)      mm_set_coordinate(typecode)

#define mm_set_real(typecode)        ((*typecode)[2]='R')
#define mm_set_pattern(typecode)     ((*typecode)[2]='P')
#define mm_set_integer(typecode)     ((*typecode)[2]='I')


#define mm_set_symmetric(typecode)  ((*typecode)[3]='S')
#define mm_set_general(typecode)    ((*typecode)[3]='G')

#define mm_clear_typecode(typecode) ((*typecode)[0]=(*typecode)[1]= \
                                     (*typecode)[2]=' ',(*typecode)[3]='G')

#define mm_initialize_typecode(typecode) mm_clear_typecode(typecode)


/********************* Matrix Market error codes ***************************/


#define MM_COULD_NOT_READ_FILE	11
#define MM_PREMATURE_EOF	12
#define MM_NOT_MTX		13
#define MM_NO_HEADER		14
#define MM_UNSUPPORTED_TYPE	15
#define MM_LINE_TOO_LONG	16
#define MM_COULD_NOT_WRITE_FILE	17


/******************** Matrix Market internal definitions ********************

   MM_matrix_typecode: 4-character sequence
        ojbect          sparse/   data        storage
        dense           type      scheme

   string position:      [0]        [1]			[2]         [3]

   Matrix typecode:  M(atrix)   C(oord)     R(eal)      G(eneral)
                     A(array)	P(attern)  S(ymmetric)  I(nteger)

 ***********************************************************************/

#define MM_MTX_STR        "matrix"
#define MM_ARRAY_STR      "array"
#define MM_DENSE_STR      "array"
#define MM_COORDINATE_STR "coordinate"
#define MM_SPARSE_STR     "coordinate"
#define MM_REAL_STR       "real"
#define MM_INT_STR        "integer"
#define MM_GENERAL_STR    "general"
#define MM_SYMM_STR       "symmetric"
#define MM_PATTERN_STR    "pattern"


/*  high level routines */
int mm_write_mtx_crd(char fname[], int M, int N, int nz, int Itab[], int Jtab[],
                     double val[], MM_typecode matcode);
int mm_read_mtx_crd_data(FILE *f, int M, int N, int nz, int Itab[], int Jtab[],
                         double val[], MM_typecode matcode);
int mm_read_mtx_crd_entry(FILE *f, int *Itab, int *Jtab, double *real, double *img,
                          MM_typecode matcode);

int mm_read_unsymmetric_sparse(const char *fname, int *M_, int *N_, int *nz_,
                double **val_, int **I_, int **J_);

#endif /* _mmio_h_ */

