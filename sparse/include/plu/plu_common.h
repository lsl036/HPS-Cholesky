#include "plu_api.h"
#include "plu_spm.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <unistd.h>
#include <errno.h>
#include <limits.h>


#if defined(PLU_OS_WINDOWS)
#include <windows.h>
#define COMMON_RANDOM_RAND 1
#endif


#ifndef _plu_config_h_
#define _plu_config_h_

#define PLU_VERSION_MAJOR 
#define PLU_VERSION_MINOR 
#define PLU_VERSION_MICRO 

/* system */
#define HAVE_PTHREAD
#define HAVE_SCHED_SETAFFINITY
#define HAVE_CLOCK_GETTIME

#define HAVE_ATOMIC_GCC_32_BUILTINS
#define HAVE_ATOMIC_GCC_64_BUILTINS

#define HAVE_STDARG_H
#define HAVE_UNISTD_H
#define HAVE_VA_COPY
#define HAVE_GETOPT_LONG
#define HAVE_GETRUSAGE
#define HAVE_GETOPT_H
#define HAVE_ERRNO_H
#define HAVE_STDDEF_H
#define HAVE_LIMITS_H
#define HAVE_STRING_H
#define HAVE_COMPLEX_H
#define HAVE_GETLINE
#define HAVE_MKDTEMP

/* Architecture */

/* Ordering options */
#define PLU_ORDERING_METIS

/* Symbolic factorization options */

/* Analyze options */
#define PLU_BLEND_DEEPEST_DISTRIB

/* #undef FORGET_PARTITION */

/* Numerical factorization options */

/* Models */

/* Datatypes used */
#define PLU_INT64

#if defined(HAVE_FALLTHROUGH)
#define plu_attr_fallthrough __attribute__((fallthrough))
#else
#define plu_attr_fallthrough do {} while(0)
#endif

#if defined(WIN32) || defined(_WIN32)
#define PLU_OS_WINDOWS 1
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

#endif /* _plu_config_h_ */

#ifndef _plu_h_
#define _plu_h_


#ifndef _plu_datatypes_h_
#define _plu_datatypes_h_

#define PLU_INT_MAX LONG_MAX


/** 
 * Double that are not converted through precision generator functions
 **/
typedef double plu_fixdbl_t;

/** 
 * Complex numbers (Extracted from PaRSEC project)
 **/
#if defined(_MSC_VER) && !defined(__INTEL_COMPILER)
/* Windows and non-Intel compiler */
#include <complex>
typedef std::complex<float>  plu_complex32_t;
typedef std::complex<double> plu_complex64_t;
#else
typedef float  _Complex      plu_complex32_t;
typedef double _Complex      plu_complex64_t;
#endif

static inline size_t
plu_size_of(plu_coeftype_t type)
{
    switch(type) {
    case PluDouble:    return   sizeof(double);
    default:
        fprintf(stderr, "plu_size_of: invalid type parameter\n");
        assert(0);
        return sizeof(double);
    }
}

/**
 * Plu data structures
 **/

/* Main structure of the plu solver associated to a given problem */
struct plu_data_s;
typedef struct plu_data_s plu_data_t;

/* Graph structure (No values) */
struct plu_graph_s;
typedef struct plu_graph_s plu_graph_t;

/* Ordering structure */
struct plu_order_s;
typedef struct plu_order_s plu_order_t;

/* Solver matrix structure to store L(U)*/
struct solver_matrix_s;
typedef struct solver_matrix_s SolverMatrix;

#endif /* _plu_datatypes_h_ */

/*
 * Main function for compatibility with former versions
 */
int plu( plu_data_t **plu_data,
            int        plu_comm,
            Sparse_long    n,
            Sparse_long   *colptr,
            Sparse_long   *row,
            void           *avals,
            Sparse_long   *perm,
            Sparse_long   *invp,
            void           *b,
            Sparse_long    nrhs,
            Sparse_long   *iparm,
            double         *dparm );

/*
 * Solver initialization
 */
void pluInitParam( Sparse_long   *iparm,
                      double         *dparm );
void pluInit     ( plu_data_t **plu_data,
                      Sparse_long   *iparm,
                      double         *dparm );
void pluInitWithAffinity( plu_data_t **plu_data,
                             Sparse_long   *iparm,
                             double         *dparm,
                             const int      *bindtab );
void pluFinalize ( plu_data_t **plu_data );

/*
 * Main steps of the solver
 */
int plu_task_analyze( plu_data_t      *plu_data,
                         sparse_csc         *spm );
int plu_task_numfact( plu_data_t      *plu_data,
                         sparse_csc         *spm );
int plu_task_solve  ( plu_data_t      *plu_data,
                         Sparse_long        nrhs,
                         void               *b,
                         Sparse_long        ldb );
int plu_task_refine( plu_data_t *plu_data,
                        Sparse_long n, Sparse_long nrhs,
                        void *b, Sparse_long ldb,
                        void *x, Sparse_long ldx );

/*
 * Analyze subtasks
 */
int plu_subtask_order     ( plu_data_t      *plu_data,
                               const sparse_csc   *spm,
                               plu_order_t     *myorder );
int plu_subtask_symbfact  ( plu_data_t      *plu_data );
int plu_subtask_reordering( plu_data_t      *plu_data );
int plu_subtask_blend     ( plu_data_t      *plu_data );

/*
 * Numerical factorization subtasks
 */
int plu_subtask_spm2bcsc  ( plu_data_t      *plu_data,
                               sparse_csc         *spm );
int plu_subtask_bcsc2ctab ( plu_data_t      *plu_data );
int plu_subtask_sopalin   ( plu_data_t      *plu_data );

/*
 * Numerical solve and refinement subtasks
 */
int plu_subtask_applyorder( plu_data_t    *plu_data,
                               plu_coeftype_t flttype,
                               plu_dir_t      dir,
                               Sparse_long      m,
                               Sparse_long      n,
                               void             *b,
                               Sparse_long      ldb );
int plu_subtask_trsm( plu_data_t    *plu_data,
                         plu_coeftype_t flttype,
                         plu_side_t     side,
                         plu_uplo_t     uplo,
                         plu_trans_t    trans,
                         plu_diag_t     diag,
                         Sparse_long      nrhs,
                         void             *b,
                         Sparse_long      ldb );
int plu_subtask_diag( plu_data_t    *plu_data,
                         plu_coeftype_t flttype,
                         Sparse_long      nrhs,
                         void             *b,
                         Sparse_long      ldb );
int plu_subtask_solve( plu_data_t *plu_data,
                          Sparse_long   nrhs,
                          void          *b,
                          Sparse_long   ldb );
int plu_subtask_refine( plu_data_t *plu_data,
                           Sparse_long n, Sparse_long nrhs,
                           const void *b, Sparse_long ldb,
                                 void *x, Sparse_long ldx );
/*
 * Schur complement manipulation routines.
 */
void pluSetSchurUnknownList( plu_data_t       *plu_data,
                                Sparse_long         n,
                                const Sparse_long  *list );
int  pluGetSchur           ( const plu_data_t *plu_data,
                                void                *S,
                                Sparse_long         lds );

void pluExpand            ( const plu_data_t *plu_data,
                               sparse_csc          *spm );

/*
 * Function to provide access to the diagonal
 */
int  pluGetDiag( const plu_data_t *plu_data,
                    void                *D,
                    Sparse_long         incD );

/*
 * Function to provide a common way to read binary options in examples/testings
 */
void pluGetOptions( int argc, char **argv, int *method,
                        Sparse_long *iparm, double *dparm,
                        int *check, char **filename );

#endif /* _plu_h_ */

#include "plu_atomic.h"

#ifndef _memory_h_
#define _memory_h_

/*
 * Function: plu_protected_malloc
 *
 * PowerPC architectures don't support malloc(0). This function
 * prints a warning when it happens to avoid segfault.
 */
static inline void *plu_malloc_func( size_t size,
                                        char *filename,
                                        int line )
{
    if (size > 0) {
	return malloc(size);
    }
    else {
	fprintf(stderr, "Pb Alloc 0 %s:%d\n", filename, line);
	return (void *)NULL;
    }
}

#if defined(PLU_ARCH_PPC)
#  define memAlloc(size) plu_malloc_func(size, __FILE__, __LINE__)
#else
#  define memAlloc(size) malloc(size)
#endif

#define memFree(ptr) free((void*)(ptr))
#define memFree_null(ptr) do			\
	{					\
	    memFree( ptr );			\
	    (ptr) = NULL;			\
	} while(0)

#define MALLOC_INTERN(ptr, size, type)                          \
    do {                                                        \
        ptr = (type*)memAlloc((size) * sizeof(type)) ;          \
    } while(0)

#define MALLOC_EXTERN(ptr, size, type)		\
    ptr = (type*)malloc((size) * sizeof(type))

#define MALLOC_ERROR( _str_ )                                           \
    {                                                                   \
        fprintf(stderr, "%s allocation (line=%d,file=%s)\n",(_str_),__LINE__,__FILE__); \
        exit(-1);                                                       \
    }

/*
 * Macro: MALLOC_INTOREXTERN
 *
 * Choose between <MALLOC_INTERN> and <MALLOC_EXTERN>
 * following flag_int.
 *
 * Parameters:
 *   ptr      - address where to allocate.
 *   size     - Number of elements to allocate.
 *   types    - Type of the elements to allocate.
 *   flag_int - 1 for internal allocation, 0 for external.
 */
#define MALLOC_INTOREXTERN(ptr, size, type, flag_int) \
  do {                                                \
    if (flag_int == 1)                          \
      {                                               \
        MALLOC_INTERN(ptr, size, type);               \
      }                                               \
    else                                              \
      {                                               \
        MALLOC_EXTERN(ptr, size, type);               \
      }                                               \
  } while (0)

#define FREE_NULL_INTOREXT(ptr, flag_int)         \
  do {                                            \
    if (flag_int == 1)                      \
      {                                           \
        memFree_null(ptr);                        \
      }                                           \
    else                                          \
      {                                           \
        free(ptr);                                \
        ptr = NULL;                               \
      }                                           \
  } while (0)

#define memRealloc realloc

#endif /* _memory_h_ */

#ifndef _integer_h_
#define _integer_h_

#ifndef MIN
#  define MIN(x,y) (((x)<(y))?(x):(y))
#endif

#ifndef MAX
#  define MAX(x,y) (((x)<(y))?(y):(x))
#endif

int          intLoad     (FILE * const, Sparse_long * const);
int          intSave     (FILE * const, const Sparse_long);
void         intAscn     (Sparse_long * restrict const, const Sparse_long, const Sparse_long);
void         intPerm     (Sparse_long * restrict const, const Sparse_long);
void         intRandInit (void);
void         intSort1asc1(void * const, const Sparse_long);
void         intSort2asc1(void * const, const Sparse_long);
void         intSort2asc2(void * const, const Sparse_long);

/*
 Function: qsortIntFloatAsc

 Sort 2 arrays simultaneously, the first array is an
 array of Sparse_long and used as key for sorting.
 The second array is an array of PLU_FLOAT.

 Parameters:
 pbase       - Array of pointers to the first element of each array to sort.
 total_elems - Number of element in each array.

 Returns:
 Nothing

 */
void d_qsortIntFloatAsc(void ** const pbase,
                        const Sparse_long     total_elems);

/*
 Function: qsort2IntFloatAsc

 Sort 3 arrays simultaneously, the first array is an
 array of Sparse_long and used as primary key for sorting.
 The second array is an other array of Sparse_long used
 as secondary key.
 The third array is an array of PLU_FLOAT.

 Parameters:
 pbase       - Array of pointers to the first element of each array to sort.
 total_elems - Number of element in each array.

 Returns:
 Nothing

 */
void d_qsort2IntFloatAsc(void ** const pbase,
                         const Sparse_long     total_elems);


/*
 Function: qsort2IntAsc

 Sort 2 arrays simultaneously, the first array is an
 array of Sparse_long and used as primary key for sorting.
 The second array is an other array of Sparse_long used
 as secondary key.

 Parameters:
 pbase       - Array of pointers to the first element of each array to sort.
 total_elems - Number of element in each array.

 Returns:
 Nothing

 */
void qsort2IntAsc(void ** const pbase,
                  const Sparse_long     total_elems);

/*
 Function: qsort3IntAsc

 Sort 3 arrays simultaneously, the first array is an
 array of Sparse_long and used as primary key for sorting.
 The second and third arrays are sorted based on the primary key.

 Parameters:
 pbase       - Array of pointers to the first element of each array to sort.
 total_elems - Number of element in each array.

 Returns:
 Nothing

 */
void qsort3IntAsc(void ** const pbase,
                  const Sparse_long     total_elems);

/*
 Function: qsort2SmallIntAsc

 Sort 2 arrays simultaneously, the first array is an
 array of integers (int) and used as primary key for sorting.
 The second array is an other array of int used
 as secondary key.

 Parameters:
 pbase       - Array of pointers to the first element of each array to sort.
 total_elems - Number of element in each array.

 Returns:
 Nothing

 */
void qsort2SmallIntAsc(void ** const pbase,
                       const Sparse_long     total_elems);


Sparse_long
plu_intset_union(       Sparse_long  n1,
                           const Sparse_long *set1,
                           Sparse_long  n2,
                           const Sparse_long *set2,
                           Sparse_long *set );

static inline Sparse_long plu_imin( Sparse_long a, Sparse_long b) {
    return ( a < b ) ? a : b;
}

static inline Sparse_long plu_imax( Sparse_long a, Sparse_long b) {
    return ( a > b ) ? a : b;
}

static inline Sparse_long plu_iceil( Sparse_long a, Sparse_long b) {
    return ( a + b - 1 ) / b;
}

Sparse_long *plu_int_convert( Sparse_long n, int *input );

#endif /* _integer_h_ */

#include "plu_symbol.h"

#ifndef _common_h_
#define _common_h_

#ifndef _timing_h_
#define _timing_h_


typedef double Clock;

/** TIMING SYSTEM-SPECIFICS **/
#if defined(HAVE_CLOCK_GETTIME)
#include <unistd.h>
#include <time.h>
static inline double clockGet(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    return ((double) ts.tv_sec + (double) ts.tv_nsec * (double)1.0e-9L);
}

#elif defined(HAVE_GETRUSAGE)
#include <sys/time.h>
static inline double clockGet(void)
{
    struct rusage       data;
    getrusage (RUSAGE_SELF, &data);
    return (((double) data.ru_utime.tv_sec  + (double) data.ru_stime.tv_sec) +
            ((double) data.ru_utime.tv_usec + (double) data.ru_stime.tv_usec) *
            1.0e-6L);
}

#elif defined(__IA64)
static inline double clockGet(void)
{
    uint64_t ret;
    __asm__ __volatile__ ("mov %0=ar.itc" : "=r"(ret));
    return (double)ret;
}

#elif defined(__X86)
static inline double clockGet(void)
{
    uint64_t ret;
    unsigned hi, lo;
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    ret = ( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 );
    return (double)ret;
}

#elif defined(__bgp__)

#include <bpcore/ppc450_inlines.h>
static inline double clockGet(void)
{
    return (double)_bgp_GetTimeBase();
}

#elif defined(HAVE_GETRUSAGE)

#include <sys/time.h>
static inline double clockGet(void)
{
    struct rusage       data;
    getrusage (RUSAGE_SELF, &data);
    return (((double) data.ru_utime.tv_sec  + (double) data.ru_stime.tv_sec) +
            ((double) data.ru_utime.tv_usec + (double) data.ru_stime.tv_usec) *
            1.0e-6L);
}

#else

#include <sys/time.h>
static inline double clockGet(void)
{
    struct timeval tv;
    gettimeofday( &tv, NULL );
    return ((double)tv.tv_sec + (double) tv.tv_usec * (double)1.0e-6L);
}

#endif

#define clockInit(clk)  do { clk = clockGet(); } while(0)
#define clockVal(clk)   (clk)

#define clockStart(clk) do { clk = clockGet(); } while(0)
#define clockStop(clk)  do { clk = clockGet() - (clk); } while(0)

#define clockSyncStart(clk) do { clk = clockGet(); } while(0)
#define clockSyncStop(clk)  do { clk = clockGet() - (clk); } while(0)

#endif /* _timing_h_ */


/*
  Macro: EXIT
  将IPARM_ERROR_NUMBER设置为module + error，转储参数并退出。
参数：
     module-发生错误的模块。
     错误-设置IPARM_ERROR_NUMBER的值。
*/
#if defined(PLU_DEBUG_EXIT_ON_SIGSEGV)
#define EXIT(module,error) { *(int *)0 = error; }
#else
#define EXIT(module,error) { abort(); }
#endif

/********************************************************************
 * CBLAS value address
 */
#ifndef CBLAS_SADDR
#define CBLAS_SADDR( a_ ) (&(a_))
#endif

/*
 * Get environment variable
 */
#if defined(PLU_OS_WINDOWS)

static inline int
plu_setenv( const char *var, const char *value, int overwrite ) {
    return !(SetEnvironmentVariable( var, value ));
}

static inline char *
plu_getenv( const char *var ) {
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
plu_cleanenv( char *str ) {
    if (str != NULL) free(str);
}

#else /* Other OS systems */

static inline int
plu_setenv( const char *var, const char *value, int overwrite ) {
    return setenv( var, value, overwrite );
}

static inline char *
plu_getenv( const char *var ) {
    return getenv( var );
}

static inline void
plu_cleanenv( char *str ) {
    (void)str;
}

#endif


static inline int
plu_env_is_set_to(char * str, char * value) {
    char * val;
    if ( (val = plu_getenv(str)) &&
         !strcmp(val, value))
        return 1;
    return 0;
}

static inline int
plu_env_is_on(char * str) {
    return plu_env_is_set_to(str, "1");
}

static inline
int plu_getenv_get_value_int(char * string, int default_value) {
    long int ret;
    int base = 10;
    char *endptr;
    char *str = plu_getenv(string);
    if (str == NULL) return default_value;

    ret = strtol(str, &endptr, base);

    /* Check for various possible errors */
    if ((errno == ERANGE && (ret == LONG_MAX || ret == LONG_MIN))
        || (errno != 0 && ret == 0)) {
        perror("strtol");
        return default_value;
    }

    if (endptr == str) {
        return default_value;
    }

    if (*endptr != '\0')        /* Not necessarily an error... */
        fprintf(stderr, "Further characters after %s value: %s\n", string, endptr);
    return (int)ret;
}

/* **************************************** */

static inline void set_iparm(Sparse_long *iparm, plu_iparm_t offset, Sparse_long value)
{
    if (iparm != NULL) iparm[offset] = (Sparse_long)value;
}

static inline void set_dparm(double *dparm, plu_dparm_t offset, double value)
{
    if (dparm != NULL) dparm[offset] = (double)value;
}

void api_dumparm(FILE *stream, Sparse_long *iparm, double *dparm);

#if !defined(HAVE_GETLINE)
ssize_t getdelim(char **buf, size_t *bufsiz, int delimiter, FILE *fp);
ssize_t getline(char **buf, size_t *bufsiz, FILE *fp);
#endif

#endif /* _common_h_ */

#ifndef _bcsc_h_
#define _bcsc_h_

/**
 * @brief Compressed colptr format for the bcsc
 */
typedef struct bcsc_cblk_s {
    Sparse_long  colnbr; /**< Number of columns in the block column.                                    */
    Sparse_long *coltab; /**< Array of indexes of the start of each column in the row and value arrays. */
} bcsc_cblk_t;

/**
 * @brief Internal column block distributed CSC matrix.
 */
struct plu_bcsc_s {
    int           gN;      /**< Global number of vertices                                                      */
    int           n;       /**< Local number of vertices                                                       */
    int           mtxtype; /**< Matrix structure: PluGeneral, PluSymmetric.           */
    int           flttype; /**< valtab datatype: PluDouble */
    Sparse_long     cscfnbr; /**< Number of column blocks.                                                       */
    bcsc_cblk_t   *cscftab; /**< Array of Block column structures of size cscfnbr. (<plu_bcscFormat_t>)      */
    Sparse_long     *rowtab;  /**< Array of rows in the matrix.                                                   */
    void          *Lvalues; /**< Array of values of the matrix A                                                */
    void          *Uvalues; /**< Array of values of the matrix A^t                                              */
};

double bcscInit( const sparse_csc     *spm,
                 const plu_order_t *ord,
                 const SolverMatrix   *solvmtx,
                 Sparse_long          initAt,
                 plu_bcsc_t        *bcsc );

void   bcscExit( plu_bcsc_t *bcsc );

Sparse_long
bcsc_init_centralized_coltab( const sparse_csc     *spm,
                              const plu_order_t *ord,
                              const SolverMatrix   *solvmtx,
                                    plu_bcsc_t  *bcsc );

void
bcsc_init_centralized( const sparse_csc     *spm,
                       const plu_order_t *ord,
                       const SolverMatrix   *solvmtx,
                             Sparse_long    initAt,
                             plu_bcsc_t  *bcsc );

void
bcsc_restore_coltab( plu_bcsc_t *bcsc );

#endif /* _bcsc_h_ */

#ifndef _bvec_h_
#define _bvec_h_

void *bvec_malloc( size_t size );

void  bvec_free( void *x );

#endif /* _bvec_h_ */


#ifndef _bcsc_d_h_
#define _bcsc_d_h_

void bcsc_dinit_centralized( const sparse_csc     *spm,
                             const plu_order_t *ord,
                             const SolverMatrix   *solvmtx,
                             const Sparse_long   *col2cblk,
                                   int             initAt,
                                   plu_bcsc_t  *bcsc );

void bvec_daxpy_seq( plu_data_t            *plu_data,
                     Sparse_long              n,
                     double        alpha,
                     const double *x,
                     double       *y );
void bvec_daxpy_smp( plu_data_t            *plu_data,
                     Sparse_long              n,
                     double        alpha,
                     const double *x,
                     double       *y );

void bvec_dcopy_seq( plu_data_t            *plu_data,
                     Sparse_long              n,
                     const double *x,
                     double       *y );
void bvec_dcopy_smp( plu_data_t            *plu_data,
                     Sparse_long              n,
                     const double *x,
                     double       *y );

#if defined(PRECISION_z) || defined(PRECISION_c)
double bvec_ddot_seq( plu_data_t            *plu_data,
                                   Sparse_long              n,
                                   const double *x,
                                   const double *y );
double bvec_ddot_smp( plu_data_t            *plu_data,
                                   Sparse_long              n,
                                   const double *x,
                                   const double *y );
#endif

double bvec_ddot_seq( plu_data_t            *plu_data,
                                   Sparse_long              n,
                                   const double *x,
                                   const double *y );
double bvec_ddot_smp( plu_data_t            *plu_data,
                                   Sparse_long              n,
                                   const double *x,
                                   const double *y );

void bvec_dgemv_seq( plu_data_t            *plu_data,
                     Sparse_long              m,
                     Sparse_long              n,
                     double        alpha,
                     const double *A,
                     Sparse_long              lda,
                     const double *x,
                     double        beta,
                     double       *y );
void bvec_dgemv_smp( plu_data_t            *plu_data,
                     Sparse_long              m,
                     Sparse_long              n,
                     double        alpha,
                     const double *A,
                     Sparse_long              lda,
                     const double *x,
                     double        beta,
                     double       *y );

double bvec_dnrm2_seq( plu_data_t            *plu_data,
                       Sparse_long              n,
                       const double *x );
double bvec_dnrm2_smp( plu_data_t            *plu_data,
                       Sparse_long              n,
                       const double *x );

void bvec_dscal_seq( plu_data_t      *plu_data,
                     Sparse_long        n,
                     double  alpha,
                     double *x );
void bvec_dscal_smp( plu_data_t      *plu_data,
                     Sparse_long        n,
                     double  alpha,
                     double *x );

int bvec_dlapmr( int thread_safe,
                 plu_dir_t        dir,
                 Sparse_long        m,
                 Sparse_long        n,
                 double *A,
                 Sparse_long        lda,
                 Sparse_long       *perm );


double bcsc_dnorm( plu_normtype_t    ntype,
                   const plu_bcsc_t *bcsc );

void bcsc_dspsv( plu_data_t      *plu_data,
                 double *b );

void bcsc_dspmv( const plu_data_t      *plu_data,
                 plu_trans_t            trans,
                 double        alpha,
                 const double *x,
                 double        beta,
                 double       *y );

void bcsc_dspmv_seq( const plu_data_t      *plu_data,
                     plu_trans_t            trans,
                     double        alpha,
                     const double *x,
                     double        beta,
                     double       *y );
void bcsc_dspmv_smp( const plu_data_t      *plu_data,
                     plu_trans_t            trans,
                     double        alpha,
                     const double *x,
                     double        beta,
                     double       *y );


#endif /* _bcsc_d_h_ */

#ifndef _out_h_
#define _out_h_

#define OUT_HEADER                                              \
    "  Schedulers:\n"                                           \
    "    sequential:                           %8s\n"           \
    "    thread static:                        %8s\n"           \
    "    thread dynamic:                       %8s\n"           \
    "  Number of threads per process:          %8d\n"           \
    "  Computational models\n"                                  \
    "    CPU: %41s\n"                                           \
    "  Low rank parameters:\n"                                  \
    "    Strategy                      %16s\n"

#define OUT_HEADER_LR                                           \
    "    Tolerance                             %8.0e\n"         \
    "    Compress method                       %8s\n"           \
    "    Compress minimal width                %8ld\n"          \
    "    Compress minimal height               %8ld\n"          \
    "    Compress min ratio                    %8f\n"           \
    "    Tolerance criterion per block         %8s\n"           \
    "    Orthogonalization method              %8s\n"           \
    "    Splitting Strategy                    %8s\n"           \
    "    Levels of projections                 %8ld\n"          \
    "    Levels of kway                        %8ld\n"          \
    "    Projections distance                  %8ld\n"          \
    "    Projections depth                     %8ld\n"          \
    "    Projections width                     %8ld\n"

#define OUT_STEP_ORDER                                          \
    "+-------------------------------------------------+\n"     \
    "  Ordering step :\n"
#define OUT_SUBSTEP_GRAPH                       \
    "    Prepare graph structure:\n"
#define OUT_ORDER_SYMGRAPH                      \
    "      Symmetrizing graph\n"
#define OUT_ORDER_NODIAG                        \
    "      Removing diagonal elements\n"
#define OUT_ORDER_SORT                          \
    "      Sort row indexes in each column\n"
#define OUT_ORDER_INIT                          \
    "    Compute ordering\n"
#define OUT_ORDER_METHOD                        \
    "    Ordering method is: %s\n"
#define OUT_ORDER_TIME                                  \
    "    Time to compute ordering              %e s\n"

#define OUT_STEP_FAX                                            \
    "+-------------------------------------------------+\n"     \
    "  Symbolic factorization step:\n"
#define OUT_FAX_METHOD                          \
    "    Symbol factorization using: %s\n"
#define OUT_FAX_SUMMARY                                                 \
    "    Number of nonzeroes in L structure    %8ld\n"                  \
    "    Fill-in of L                          %8lf\n"                  \
    "    Time to compute symbol matrix         %e s\n"


#define OUT_STEP_REORDER                                        \
    "+-------------------------------------------------+\n"     \
    "  Reordering step:\n"                                      \
    "    Split level                           %8ld\n"          \
    "    Stoping criterion                     %8ld\n"
#define OUT_REORDERING_TIME                             \
    "    Time for reordering                   %e s\n"
#define OUT_REORDERING_OPS                                              \
    "    Iops for the last supernode           %8ld ( %5.2lf%% )\n"     \
    "    Iops for the reordering               %8ld\n"

#define OUT_STEP_BLEND                                          \
    "+-------------------------------------------------+\n"     \
    "  Analyse step:\n"
#define OUT_BLEND_CONF                                  \
    "    Number of cluster                     %8ld\n"  \
    "    Number of processor per cluster       %8ld\n"  \
    "    Number of thread per MPI process      %8ld\n"

#define OUT_BLEND_CHKSMBMTX                     \
    "    Check the symbol matrix\n"
#define OUT_BLEND_CHKSOLVER                     \
    "    Check the solver matrix\n"
#define OUT_BLEND_ELIMTREE                      \
    "    Building elimination tree\n"
#define OUT_BLEND_ELIMTREE_TIME                         \
    "    Elimination tree built in             %e s\n"
#define OUT_BLEND_COSTMATRIX                    \
    "    Building cost matrix\n"
#define OUT_BLEND_COSTMATRIX_TIME                       \
    "    Cost matrix built in                  %e s\n"
#define OUT_BLEND_ELIMTREE_TOTAL_COST                   \
    "    Total estimated cost of the etree     %e s\n"
#define OUT_BLEND_PROPMAP                       \
    "    Perform proportional mapping\n"
#define OUT_BLEND_PROPMAP_TIME                          \
    "    Proportional mapping done in          %e s\n"
#define OUT_BLEND_SPLITSYMB                     \
    "    Split large symbolic blocks\n"
#define OUT_BLEND_SPLITSYMB_TIME                        \
    "    Symbol split done in                  %e s\n"
#define OUT_BLEND_BUILDSIMU                     \
    "    Build simulation structures\n"
#define OUT_BLEND_BUILDSIMU_TIME                        \
    "    Simulation structures built in        %e s\n"  \
    "    Number of tasks found                 %8ld\n"
#define OUT_BLEND_SIMU                                  \
    "    Start simulation (Data distribution)\n"
#define OUT_BLEND_SIMU_TIME                             \
    "    Simulation done in                    %e s\n"
#define OUT_BLEND_ELIMGRAPH                     \
    "    Building elimination graph\n"
#define OUT_BLEND_ELIMGRAPH_TIME                        \
    "    Elimination graph built in            %e s\n"
#define OUT_BLEND_SOLVER                        \
    "    Building solver structure\n"
#define OUT_BLEND_SOLVER_TIME                           \
    "    Solver built in                       %e s\n"
#define OUT_BLEND_TIME                                  \
    "    Time for analyze                      %e s\n"

#define OUT_BLEND_SUMMARY                                               \
    "    Number of non-zeroes in blocked L     %8ld\n"                  \
    "    Fill-in                               %8lf\n"                  \
    "    Number of operations in full-rank: %-5s    %5.2lf %cFlops\n"   \
    "    Prediction:\n"                                                 \
    "      Model                       %20s\n"                          \
    "      Time to factorize                   %e s\n"                  \
    "    Time for analyze                      %e s\n"

#define OUT_STEP_SOPALIN                                          \
    "+-------------------------------------------------+\n"     \
    "  Factorization step:\n"                                   \
    "    Factorization used: %s\n"

#define OUT_BCSC_TIME                                   \
    "    Time to initialize internal csc       %e s\n"

#define OUT_COEFTAB_TIME                                \
    "    Time to initialize coeftab            %e s\n"

#define OUT_SOPALIN_TIME                                                \
    "    Time to factorize                     %e s (%5.2lf %cFlop/s)\n" \
    "    Number of operations                       %5.2lf %cFlops\n"   \
    "    Number of static pivots               %8ld\n"

#define OUT_LOWRANK_SUMMARY                                     \
    "    Compression:\n"                                        \
    "      Elements removed             %8ld / %8ld\n"          \
    "      Memory saved              %.3g %co / %.3g %co\n"

#define OUT_MATRIX_SIZE       "  Matrix size                                   %ld x %ld\n"
#define OUT_NNZ               "  Number of nonzeros in A                       %ld\n"

#define OUT_GLOBAL_NNZL       "   Number of nonzeroes in L structure      %ld\n"
#define OUT_GLOBAL_FILLIN     "   Fill-in                                 %lf\n"
#define OUT_GLOBAL_THFLOPCNT  "   Number of theoretical flop            %.5g %cflops\n"
#define OUT_GLOBAL_RLFLOPCNT  "   Number of performed flop              %.5g %cflops\n"

#define TIME_TO_ANALYSE       "   Time to analyze                              %.3g s\n"
#define NNZERO_WITH_FILLIN_TH "   Number of nonzeros in factorized matrix      %ld\n"
#define NNZERO_WITH_FILLIN    "%d : Number of nonzeros (local block structure) %ld\n"
#define SOLVMTX_WITHOUT_CO    "%d : SolverMatrix size (without coefficients)   %.3g %s\n"
#define OUT_FILLIN_TH         "   Fill-in                                      %lg\n"
#define NUMBER_OP_LU          "   Number of operations (LU)                    %g\n"
#define TIME_FACT_PRED        "   Prediction Time to factorize (%s) %.3g s\n"
#define OUT_COEFSIZE          "   Maximum coeftab size (cefficients)           %.3g %co\n"
#define OUT_REDIS_CSC         "   Redistributing user CSC into PLU distribution\n"
#define OUT_REDIS_RHS         "   Redistributing user RHS into PLU distribution\n"
#define OUT_REDIS_SOL         "   Redistributing solution into Users' distribution\n"
#define OUT2_SOP_BINITG       "   --- Sopalin : Allocation de la structure globale ---\n"
#define OUT2_SOP_EINITG       "   --- Fin Sopalin Init                             ---\n"
#define OUT2_SOP_TABG         "   --- Initialisation des tableaux globaux          ---\n"
#define OUT2_SOP_BINITL       "   --- Sopalin : Local structure allocation         ---\n"
#define OUT2_SOP_NOTBIND      "   --- Sopalin : Threads are NOT binded             ---\n"
#define OUT2_SOP_BIND         "   --- Sopalin : Threads are binded                 ---\n"
#define OUT2_FUN_STATS        "     - %3ld : Envois %5ld - Receptions %5ld          -\n"
#define OUT2_SOP_BSOP         "   --- Sopalin Begin                                ---\n"
#define OUT2_SOP_ESOP         "   --- Sopalin End                                  ---\n"
#define OUT4_UPDO_TIME_INIT   " [%d][%d] Solve initialization time : %lg s\n"
#define OUT4_UPDO_COMM_TIME   " [%d][%d] Solve communication time : %lg s\n"
#define OUT4_FACT_COMM_TIME   " [%d][%d] Factorization communication time : %lg s\n"
#define OUT2_SOP_DOWN         "   --- Down Step                                    ---\n"
#define OUT2_SOP_DIAG         "   --- Diag Step                                    ---\n"
#define OUT2_SOP_UP           "   --- Up Step                                      ---\n"
#define GEN_RHS_1             "   Generate RHS for X=1\n"
#define GEN_RHS_I             "   Generate RHS for X=i\n"
#define GEN_SOL_0             "   Generate X0=0\n"
#define OUT_ITERREFINE_GMRES    "   GMRES :\n"
#define OUT_ITERREFINE_PIVOT    "   Simple refinement :\n"
#define OUT_ITERREFINE_BICGSTAB "   BICGSTAB :\n"
#define OUT_ITERREFINE_GRAD     "   Conjuguate gradient :\n"
#define OUT_ITERREFINE_ITER     "    - iteration %d :\n"
#define OUT_ITERREFINE_TTS      "         time to solve                          %.3g s\n"
#define OUT_ITERREFINE_TTT      "         total iteration time                   %.3g s\n"
#define OUT_ITERREFINE_ERR      "         error                                  %.5g\n"
#define OUT_ITERREFINE_NORMA    "         ||A||                                  %.5g\n"
#define OUT_ITERREFINE_NORMR    "         ||r||                                  %.5g\n"
#define OUT_ITERREFINE_NORMB    "         ||b||                                  %.5g\n"
#define OUT_ITERREFINE_BDIVR    "         ||r||/||b||                            %.5g\n"
#define OUT_REDISCSCDTIME     "   Time to redistribute cscd                    %.3g s\n"
#define OUT_FILLCSCTIME       "   Time to fill internal csc                    %.3g s\n"
#define OUT_MAX_MEM_AF_SOP    "   Max memory used after factorization          %.3g %s\n"
#define OUT_MEM_USED_AF_SOP   "   Memory used after factorization              %.3g %s\n"
#define MAX_MEM_AF_CL         "   Max memory used after clean                  %.3g %s\n"
#define MEM_USED_AF_CL        "   Memory used after clean                      %.3g %s\n"
#define OUT_STATIC_PIVOTING   "   Static pivoting                              %ld\n"
#define OUT_ESP_NBTASKS       "   Number of tasks added by esp                 %ld\n"
#define OUT_TIME_FACT         "   Time to factorize                            %.3g s  (%.3g %s)\n"
#define OUT_FLOPS_FACT        "   FLOPS during factorization                   %.5g %s\n"
#define OUT_TIME_SOLV         "    Time to solve                         %e s\n"
#define OUT_REFINE_ITER_NORM  "    Refinement                            %ld iterations, norm=%e\n"
#define OUT_PREC1             "    ||b-Ax||/||b||                        %e\n"
#define OUT_PREC2             "    max_i(|b-Ax|_i/(|b| + |A||x|)_i       %e\n"
#define OUT_TIME_REFINE       "    Time for refinement                   %e s\n"
#define OUT_END               " +--------------------------------------------------------------------+\n"

/*
 * Printing function to redirect to the correct file
 */
#if defined(__GNUC__)
static inline void plu_print(const char *fmt, ...) __attribute__((format(printf,1,2)));
static inline void plu_print_error  ( const char *fmt, ...) __attribute__((format(printf,1,2)));
static inline void plu_print_warning( const char *fmt, ...) __attribute__((format(printf,1,2)));
#endif

static inline void
plu_print(const char *fmt, ...)
{
    // va_list ap;
    // va_start(ap, fmt);
    // vfprintf(stdout, fmt, ap );
    // va_end(ap);
}

static inline void
plu_print_error( const char *fmt, ... )
{
    va_list arglist;
    va_start(arglist, fmt);
    vfprintf(stderr, fmt, arglist);
    va_end(arglist);
    // assert(0);
}

static inline void
plu_print_warning( const char *fmt, ... )
{
    va_list arglist;
    va_start(arglist, fmt);
    fprintf(stderr, "WARNING: ");
    vfprintf(stderr, fmt, arglist);
    va_end(arglist);
}

#define errorPrint  plu_print_error
#define errorPrintW plu_print_warning

static inline double
plu_print_value( double flops )
{
    static double ratio = (double)(1<<10);
    int unit = 0;

    while ( (flops > ratio) && (unit < 9) ) {
        flops /= ratio;
        unit++;
    }
    return flops;
}

static inline char
plu_print_unit( double flops )
{
    static char units[9] = { ' ', 'K', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y' };
    static double ratio = (double)(1<<10);
    int unit = 0;

    while ( (flops > ratio) && (unit < 9) ) {
        flops /= ratio;
        unit++;
    }
    return units[unit];
}

static inline char *
pluFactotypeStr( plu_factotype_t ft ) {
    switch( ft ) {
    case PluFactLU:
        return "LU";
     default:
        return "None";
    }
}

void   plu_gendirectories( plu_data_t *plu_data );
FILE * plu_fopenw( const char *dirname,
                      const char *filename,
                      const char *mode );
FILE * plu_fopen ( const char *filename );

#endif /* _out_h_ */

