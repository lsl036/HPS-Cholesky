
#ifndef _plu_lowrank_h_
#define _plu_lowrank_h_

/**
 * @brief List of short names for the compression kernels
 */
extern const char *compmeth_shnames[PluCompressMethodNbr];

/**
 * @brief List of long names for the compression kernels
 */
extern const char *compmeth_lgnames[PluCompressMethodNbr];

/**
 * @brief Macro to specify if the U part of a low-rank matrix is orthogonal or not (Used in LRMM functions).
 */
#define PLU_LRM3_ORTHOU (1 << 0)
/**
 * @brief Macro to specify if the U part of a low-rank matrix has been allocated and need to be freed or not (Used in LRMM functions).
 */
#define PLU_LRM3_ALLOCU (1 << 1)
/**
 * @brief Macro to specify if the V part of a low-rank matrix has been allocated and need to be freed or not (Used in LRMM functions).
 */
#define PLU_LRM3_ALLOCV (1 << 2)
/**
 * @brief Macro to specify if the the operator on B, still needs to be applied to the V part of the low-rank matrix or not (Used in LRMM functions).
 */
#define PLU_LRM3_TRANSB (1 << 3)

/**
 * @brief 定义我们接受的将矩阵压缩为低秩形式或不压缩的最小比率。
 */
extern double plu_lr_minratio;

/**
 * @brief 定义正交化方法。
 */
extern Sparse_long plu_lr_ortho;

/**
 * @brief 计算用于测试的给定矩阵大小的最大秩
 * @param[in] M The number of rows of the matrix
 * @param[in] N The number of columns of the matrix
 * @return The maximal rank accepted for this matrix size.
 */
static inline Sparse_long
core_get_rklimit_max( Sparse_long M, Sparse_long N ) {
    return plu_imin( M, N );
}

/**
 * @brief 计算在即时生产策略中给定矩阵大小下可接受的最大秩
 * @param[in] M The number of rows of the matrix
 * @param[in] N The number of columns of the matrix
 * @return The maximal rank accepted for this matrix size.
 */
static inline Sparse_long
core_get_rklimit_end( Sparse_long M, Sparse_long N ) {
    return ( plu_lr_minratio * plu_imin( M, N ) ) / 4;
}

/**
 * @brief 计算最小存储策略中给定矩阵大小可接受的最大秩
 * @param[in] M The number of rows of the matrix
 * @param[in] N The number of columns of the matrix
 * @return The maximal rank accepted for this matrix size.
 */
static inline Sparse_long
core_get_rklimit_begin( Sparse_long M, Sparse_long N ) {
    return ( plu_lr_minratio * M * N ) / ( M + N );
}

static inline Sparse_long
core_get_rklimit_test( Sparse_long M, Sparse_long N ) {
    return plu_imin( M, N );
}

extern Sparse_long (*core_get_rklimit)( Sparse_long, Sparse_long );

struct plu_lr_s;
typedef struct plu_lr_s plu_lr_t;

/**
 * @brief 以低秩形式保存矩阵的块低秩结构
 */
typedef struct plu_lrblock_s {
    int   rk;    /**< Rank of the low-rank matrix: -1 is dense, otherwise rank-rk matrix           */
    int   rkmax; /**< Leading dimension of the matrix u                                            */
    void *u;     /**< Contains the dense matrix if rk=-1, or the u factor from u vT representation */
    void *v;     /**< Not referenced if rk=-1, otherwise, the v factor                             */
} plu_lrblock_t;

/**
 * @brief 将稠密块压缩为低秩形式的函数类型。
 */
typedef plu_fixdbl_t (*fct_ge2lr_t)( int, plu_fixdbl_t, Sparse_long, Sparse_long, Sparse_long,
                                        const void *, Sparse_long, plu_lrblock_t * );

/**
 * @brief 指向ge2lr的多个算术和算法变体的指针数组
 */
extern const fct_ge2lr_t ge2lrMethods[PluCompressMethodNbr][4];

/**
 * @brief 将两个低秩块相加的函数类型。
 */
typedef plu_fixdbl_t (*fct_rradd_t)( const plu_lr_t *, plu_trans_t, const void *,
                                        Sparse_long, Sparse_long, const plu_lrblock_t *,
                                        Sparse_long, Sparse_long,       plu_lrblock_t *,
                                        Sparse_long, Sparse_long );

/**
 * @brief 指向rradd的多个算术和算法变量的指针数组
 */
extern const fct_rradd_t rraddMethods[PluCompressMethodNbr][4];

/**
 * @brief 结构来定义用于低秩内核及其参数的函数类型。
 */
typedef struct plu_lr_s {
    Sparse_long compress_when;       /**< When to compress in the full solver               */
    Sparse_long compress_method;     /**< Compression method                                */
    Sparse_long compress_min_width;  /**< Minimum width to compress a supernode             */
    Sparse_long compress_min_height; /**< Minimum height to compress an off-diagonal block  */
    int          use_reltol;          /**< Enable/disable relative tolerance vs absolute one */
    double       tolerance;           /**< Absolute compression tolerance                    */
    fct_rradd_t  core_rradd;          /**< Recompression function                            */
    fct_ge2lr_t  core_ge2lr;          /**< Compression function                              */
} plu_lr_t;

/**
 * @brief 枚举用于定义块的类型。
 */
typedef enum memory_stats_e {
    FR_InDiag  = 0, /**< Full-rank block, inside a diagonal block from non-split partition*/
    FR_OffDiag = 1, /**< Full-rank block, outside a diagonal block from non-split partition*/
    LR_InDiag  = 2, /**< Low-rank block, inside a diagonal block from non-split partition*/
    LR_OffDiag = 3, /**< Low-rank block, outside a diagonal block from non-split partition*/
    LR_InSele  = 4, /**< Non-compressible block, outside a diagonal block from non-split partition*/
    LR_OffSele = 5, /**< Non-compressible block, outside a diagonal block from non-split partition*/
    LR_ToSele  = 6  /**< Non-compressible block, inside a diagonal block from non-split partition, that contributes to a selected cblk*/
} memory_stats_t;

#endif /* _plu_lowrank_h_ */

#ifndef _plu_dlrcores_h_
#define _plu_dlrcores_h_

/**
 *    这个模块包含所有使用plu_lr_t矩阵表示的低秩内核。
 */
void core_dlralloc( Sparse_long M, Sparse_long N, Sparse_long rkmax, plu_lrblock_t *A );
void core_dlrfree ( plu_lrblock_t *A );
int  core_dlrsze  ( int copy, Sparse_long M, Sparse_long N, plu_lrblock_t *A, Sparse_long newrk, Sparse_long newrkmax, Sparse_long rklimit );
int  core_dlr2ge  ( plu_trans_t trans, Sparse_long M, Sparse_long N, const plu_lrblock_t *Alr, double *A, Sparse_long lda );

void core_dlrcpy  ( const plu_lr_t *lowrank,
                    plu_trans_t transA, double alpha,
                    Sparse_long M1, Sparse_long N1, const plu_lrblock_t *A,
                    Sparse_long M2, Sparse_long N2,       plu_lrblock_t *B,
                    Sparse_long offx, Sparse_long offy );

void core_dlrconcatenate_u( double alpha,
                            Sparse_long M1, Sparse_long N1, const plu_lrblock_t *A,
                            Sparse_long M2,                        plu_lrblock_t *B,
                            Sparse_long offx,
                            double *u1u2 );
void core_dlrconcatenate_v( plu_trans_t transA1, double alpha,
                            Sparse_long M1, Sparse_long N1, const plu_lrblock_t *A,
                                             Sparse_long N2,       plu_lrblock_t *B,
                            Sparse_long offy,
                            double *v1v2 );

double core_dlrnrm( plu_normtype_t ntype, int transV,
                    Sparse_long M, Sparse_long N,
                    const plu_lrblock_t *A );

/**
 * @brief 结构存储core_dlrmm家族函数的所有参数
 */
typedef struct core_dlrmm_s {
    const plu_lr_t      *lowrank;     /**< The lowrank structure                                                 */
    plu_trans_t          transA;      /**< Specify op(A) and is equal to PluNoTrans, PluTrans, or PluTrans */
    plu_trans_t          transB;      /**< Specify op(B) and is equal to PluNoTrans, PluTrans, or PluTrans */
    Sparse_long            M;           /**< Number of rows     of the A matrix                                    */
    Sparse_long            N;           /**< Number of columns  of the B matrix                                    */
    Sparse_long            K;           /**< Number of columns  of the A matrix (= number of rows of the B matrix) */
    Sparse_long            Cm;          /**< Number of rows     of the C matrix that receives the AB contribution  */
    Sparse_long            Cn;          /**< Number of columns  of the C matrix that receives the AB contribution  */
    Sparse_long            offx;        /**< Horizontal offsets of the AB product in the C matrix                  */
    Sparse_long            offy;        /**< Vertical   offsets of the AB product in the C matrix                  */
    double      alpha;       /**< The alpha factor                                                      */
    const plu_lrblock_t *A;           /**< The A matrix described in a low-rank structure                        */
    const plu_lrblock_t *B;           /**< The B matrix described in a low-rank structure                        */
    double      beta;        /**< The beta factor                                                       */
    plu_lrblock_t       *C;           /**< The C matrix described in a low-rank structure                        */
    double     *work;        /**< The pointer to an available workspace                                 */
    Sparse_long            lwork;       /**< The size of the given workspace                                       */
    Sparse_long            lwused;      /**< The size of the workspace that is already used                        */
    plu_atomic_lock_t   *lock;        /**< The lock to protect the concurrent accesses on the C matrix           */
} core_dlrmm_t;

/**
 * @brief 初始化core_dlrmm系列函数的所有参数，以简化访问
 */
#define PASTE_CORE_DLRMM_PARAMS(_a_)                   \
    const plu_lr_t      *lowrank = (_a_)->lowrank;  \
    plu_trans_t          transA  = (_a_)->transA;   \
    plu_trans_t          transB  = (_a_)->transB;   \
    Sparse_long            M       = (_a_)->M;        \
    Sparse_long            N       = (_a_)->N;        \
    Sparse_long            K       = (_a_)->K;        \
    Sparse_long            Cm      = (_a_)->Cm;       \
    Sparse_long            Cn      = (_a_)->Cn;       \
    Sparse_long            offx    = (_a_)->offx;     \
    Sparse_long            offy    = (_a_)->offy;     \
    double      alpha   = (_a_)->alpha;    \
    const plu_lrblock_t *A       = (_a_)->A;        \
    const plu_lrblock_t *B       = (_a_)->B;        \
    double      beta    = (_a_)->beta;     \
    plu_lrblock_t       *C       = (_a_)->C;        \
    double     *work    = (_a_)->work;     \
    Sparse_long            lwork   = (_a_)->lwork;    \
    plu_atomic_lock_t   *lock    = (_a_)->lock;

/**
 * @brief 将core_dlrmm家族函数的所有参数Void为静默警告
 */
#define PASTE_CORE_DLRMM_VOID                   \
    (void)lowrank;                              \
    (void)transA;                               \
    (void)transB;                               \
    (void)M;                                    \
    (void)N;                                    \
    (void)K;                                    \
    (void)Cm;                                   \
    (void)Cn;                                   \
    (void)offx;                                 \
    (void)offy;                                 \
    (void)alpha;                                \
    (void)A;                                    \
    (void)B;                                    \
    (void)beta;                                 \
    (void)C;                                    \
    (void)work;                                 \
    (void)lwork;                                \
    (void)lock

/**
 * @brief 如果提供的工作区指针中有可用空间，则函数获得工作区指针
 * @param[inout] params  The parameters structure for core_dlrmm family functions
 * @param[in]    newsize The required workspace size in number of elements
 * @return The pointer to the workspace if enough space available, NULL otherwise.
 */
static inline double *
core_dlrmm_getws( core_dlrmm_t *params,
                  ssize_t newsize )
{
    double *work = NULL;
    if ( (params->lwused + newsize) <= params->lwork )
    {
        work = params->work + params->lwused;
        params->lwused += newsize;
    }
    /* else */
    /* { */
    /*     if ( (params->work == NULL) || (params->lwused == 0) ) */
    /*     { */
    /*         params->work = realloc( params->work, newsize * sizeof(double) ); */
    /*         params->lwork  = newsize; */
    /*         params->lwused = newsize; */
    /*         work = params->work; */
    /*     } */
    /* } */
    return work;
}

/**
 *      @name update_fr函数，对全秩矩阵执行更新
 */
plu_fixdbl_t core_dfrfr2fr( core_dlrmm_t *params );
plu_fixdbl_t core_dfrlr2fr( core_dlrmm_t *params );
plu_fixdbl_t core_dlrfr2fr( core_dlrmm_t *params );
plu_fixdbl_t core_dlrlr2fr( core_dlrmm_t *params );

/**
 *  update_lr函数为低秩矩阵的更新准备AB乘积
 */
plu_fixdbl_t core_dfrfr2lr( core_dlrmm_t     *params,
                               plu_lrblock_t *AB,
                               int              *infomask,
                               Sparse_long      Kmax );
plu_fixdbl_t core_dfrlr2lr( core_dlrmm_t     *params,
                               plu_lrblock_t *AB,
                               int              *infomask,
                               Sparse_long      Brkmin );
plu_fixdbl_t core_dlrfr2lr( core_dlrmm_t     *params,
                               plu_lrblock_t *AB,
                               int              *infomask,
                               Sparse_long      Arkmin );
plu_fixdbl_t core_dlrlr2lr( core_dlrmm_t     *params,
                               plu_lrblock_t *AB,
                               int              *infomask );

/**
 *   add_lr函数将AB贡献以低秩格式添加到任何C矩阵
 */
plu_fixdbl_t core_dlradd( core_dlrmm_t           *params,
                             const plu_lrblock_t *AB,
                             plu_trans_t          transV,
                             int                     infomask );

plu_fixdbl_t core_dlrmm( core_dlrmm_t *params );

/**
 *    这是基于LAPACK GESVD函数的低秩核的SVD实现。
 */

plu_fixdbl_t core_dge2lr_svd( int use_reltol, plu_fixdbl_t tol, Sparse_long rklimit,
                                 Sparse_long m, Sparse_long n,
                                 const void *Avoid, Sparse_long lda, plu_lrblock_t *Alr );
plu_fixdbl_t core_drradd_svd( const plu_lr_t *lowrank, plu_trans_t transA1, const void *alphaptr,
                                 Sparse_long M1, Sparse_long N1, const plu_lrblock_t *A,
                                 Sparse_long M2, Sparse_long N2,       plu_lrblock_t *B,
                                 Sparse_long offx, Sparse_long offy);

/**
 *    这些是rank-revealing QR实现，以生成全秩矩阵的低秩表示。
 */

typedef int (*core_drrqr_cp_t)( double tol, Sparse_long maxrank, int refine, Sparse_long nb,
                                Sparse_long m, Sparse_long n,
                                double *A, Sparse_long lda,
                                Sparse_long *jpvt, double *tau,
                                double *work, Sparse_long lwork,  double *rwork );

typedef int (*core_drrqr_rt_t)( double tol, Sparse_long maxrank, Sparse_long nb,
                                Sparse_long m, Sparse_long n,
                                double *A, Sparse_long lda, double *tau,
                                double *B, Sparse_long ldb, double *tau_b,
                                double *work, Sparse_long lwork,  double normA );

plu_fixdbl_t core_dge2lr_pqrcp( int use_reltol, plu_fixdbl_t tol, Sparse_long rklimit,
                                   Sparse_long m, Sparse_long n,
                                   const void *Avoid, Sparse_long lda, plu_lrblock_t *Alr );
plu_fixdbl_t core_drradd_pqrcp( const plu_lr_t *lowrank, plu_trans_t transA1, const void *alphaptr,
                                   Sparse_long M1, Sparse_long N1, const plu_lrblock_t *A,
                                   Sparse_long M2, Sparse_long N2,       plu_lrblock_t *B,
                                   Sparse_long offx, Sparse_long offy );

plu_fixdbl_t core_dge2lr_rqrcp( int use_reltol, plu_fixdbl_t tol, Sparse_long rklimit,
                                   Sparse_long m, Sparse_long n,
                                   const void *Avoid, Sparse_long lda, plu_lrblock_t *Alr );
plu_fixdbl_t core_drradd_rqrcp( const plu_lr_t *lowrank, plu_trans_t transA1, const void *alphaptr,
                                   Sparse_long M1, Sparse_long N1, const plu_lrblock_t *A,
                                   Sparse_long M2, Sparse_long N2,       plu_lrblock_t *B,
                                   Sparse_long offx, Sparse_long offy );

plu_fixdbl_t core_dge2lr_tqrcp( int use_reltol, plu_fixdbl_t tol, Sparse_long rklimit,
                                   Sparse_long m, Sparse_long n,
                                   const void *Avoid, Sparse_long lda, plu_lrblock_t *Alr );
plu_fixdbl_t core_drradd_tqrcp( const plu_lr_t *lowrank, plu_trans_t transA1, const void *alphaptr,
                                   Sparse_long M1, Sparse_long N1, const plu_lrblock_t *A,
                                   Sparse_long M2, Sparse_long N2,       plu_lrblock_t *B,
                                   Sparse_long offx, Sparse_long offy );

plu_fixdbl_t core_dge2lr_rqrrt( int use_reltol, plu_fixdbl_t tol, Sparse_long rklimit,
                                   Sparse_long m, Sparse_long n,
                                   const void *Avoid, Sparse_long lda, plu_lrblock_t *Alr );


plu_fixdbl_t core_dge2lr_qrcp( core_drrqr_cp_t rrqrfct,
                                  int use_reltol, plu_fixdbl_t tol, Sparse_long rklimit,
                                  Sparse_long m, Sparse_long n,
                                  const void *Avoid, Sparse_long lda,
                                  plu_lrblock_t *Alr );
plu_fixdbl_t core_dge2lr_qrrt( core_drrqr_rt_t rrqrfct,
                                  int use_reltol, plu_fixdbl_t tol, Sparse_long rklimit,
                                  Sparse_long m, Sparse_long n,
                                  const void *Avoid, Sparse_long lda,
                                  plu_lrblock_t *Alr);

plu_fixdbl_t core_drradd_qr( core_drrqr_cp_t rrqrfct,
                                const plu_lr_t *lowrank, plu_trans_t transA1, const void *alphaptr,
                                Sparse_long M1, Sparse_long N1, const plu_lrblock_t *A,
                                Sparse_long M2, Sparse_long N2,       plu_lrblock_t *B,
                                Sparse_long offx, Sparse_long offy );

/**
 *  the debug routines for the low rank kernels
 */

void core_dlrdbg_printsvd( Sparse_long              M,
                           Sparse_long              N,
                           const double *A,
                           Sparse_long              lda );

int  core_dlrdbg_check_orthogonality( Sparse_long              M,
                                      Sparse_long              N,
                                      const double *A,
                                      Sparse_long              lda );

int  core_dlrdbg_check_orthogonality_AB( Sparse_long M, Sparse_long NA, Sparse_long NB,
                                         const double *A, Sparse_long lda,
                                         const double *B, Sparse_long ldb );

#endif /* _plu_dlrcores_h_ */
