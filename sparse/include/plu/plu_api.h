#ifndef _plu_api_h_
#define _plu_api_h_

/**
 * @brief Integer parameters
 */
typedef enum plu_iparm_e {
    IPARM_VERBOSE,               /**< Verbose mode (@see plu_verbose_t)                           Default: PluVerboseNo           IN  */
    IPARM_IO_STRATEGY,           /**< IO strategy  (@see plu_io_t)                                Default: PluIONo                IN  */

    /* Stats */
    IPARM_NNZEROS,               /**< Number of nonzero entries in the factorized matrix             Default: -                         OUT */
    IPARM_NNZEROS_BLOCK_LOCAL,   /**< Number of nonzero entries in the local block factorized matrix Default: -                         OUT */
    IPARM_ALLOCATED_TERMS,       /**< Maximum memory allocated for matrix terms                      Default: -                         OUT */
    IPARM_PRODUCE_STATS,         /**< Compute some statistiques (such as precision error)            Default: 0                         IN  */

    /* Scaling */
    IPARM_MC64,                  /**< MC64 operation                                                 Default: 0                         IN  */

    /* Ordering */
    IPARM_ORDERING,              /**< Choose ordering                                                Default: PluOrderMetis          IN  */
    IPARM_ORDERING_DEFAULT,      /**< Use default ordering parameters with Metis                     Default: 1                         IN  */

    /* Subset for Metis */
    IPARM_METIS_CTYPE,           /**< Metis parameters (see Metis Manual)                            Default: METIS_CTYPE_SHEM          IN  */
    IPARM_METIS_RTYPE,           /**< Metis parameters (see Metis Manual)                            Default: METIS_RTYPE_SEP1SIDED     IN  */
    IPARM_METIS_NO2HOP,          /**< Metis parameters (see Metis Manual)                            Default: 0                         IN  */
    IPARM_METIS_NSEPS,           /**< Metis parameters (see Metis Manual)                            Default: 1                         IN  */
    IPARM_METIS_NITER,           /**< Metis parameters (see Metis Manual)                            Default: 10                        IN  */
    IPARM_METIS_UFACTOR,         /**< Metis parameters (see Metis Manual)                            Default: 200                       IN  */
    IPARM_METIS_COMPRESS,        /**< Metis parameters (see Metis Manual)                            Default: 1                         IN  */
    IPARM_METIS_CCORDER,         /**< Metis parameters (see Metis Manual)                            Default: 0                         IN  */
    IPARM_METIS_PFACTOR,         /**< Metis parameters (see Metis Manual)                            Default: 0                         IN  */
    IPARM_METIS_SEED,            /**< Metis parameters (see Metis Manual)                            Default: 3452                      IN  */
    IPARM_METIS_DBGLVL,          /**< Metis parameters (see Metis Manual)                            Default: 0                         IN  */

    /* Symbolic Factorization */
    IPARM_AMALGAMATION_LVLBLAS,  /**< Amalgamation level                                             Default: 5                         IN  */
    IPARM_AMALGAMATION_LVLCBLK,  /**< Amalgamation level                                             Default: 5                         IN  */

    /* Reordering */
    IPARM_REORDERING_SPLIT,      /**< Reordering split level                                         Default: 0                         IN  */
    IPARM_REORDERING_STOP,       /**< Reordering stop criterion                                      Default: PLU_INT_MAX            IN  */

    IPARM_SPLITTING_STRATEGY,    /**< Strategy used to split supernodes                              Default: PLU_INT_MAX           IN    */
    IPARM_SPLITTING_LEVELS_PROJECTIONS,   /**< Levels of projections                                 Default: PLU_INT_MAX           IN    */
    IPARM_SPLITTING_LEVELS_KWAY ,         /**< Levels of kway                                        Default: PLU_INT_MAX           IN    */
    IPARM_SPLITTING_PROJECTIONS_DEPTH,    /**< Number of level used for projections                  Default: PLU_INT_MAX           IN    */
    IPARM_SPLITTING_PROJECTIONS_DISTANCE, /**< Distance used for projections                         Default: PLU_INT_MAX           IN    */
    IPARM_SPLITTING_PROJECTIONS_WIDTH,    /**< Width used for projections
                                           Default: PLU_INT_MAX IN */

    /* Analyze */
    IPARM_MIN_BLOCKSIZE,         /**< Minimum block size                                             Default: 160                       IN  */
    IPARM_MAX_BLOCKSIZE,         /**< Maximum block size                                             Default: 320                       IN  */
    IPARM_TASKS2D_LEVEL,         /**< 2D Distribution level (-1 for autolevel, 0 for 1D)             Default: -1                        IN  */
    IPARM_TASKS2D_WIDTH,         /**< Minimal width for 2D tasks with autolevel                      Default: IPARM_MIN_BLOCKSIZE       IN  */
    IPARM_ALLCAND,               /**< Allow all threads to be candidate in the proportional mapping  Default: 0                         IN  */

    /* Incomplete */
    IPARM_INCOMPLETE,            /**< Incomplete factorization                                       Default: 0                         IN  */
    IPARM_LEVEL_OF_FILL,         /**< Level of fill for incomplete factorization                     Default: 0                         IN  */

    /* Factorization */
    IPARM_FACTORIZATION,         /**< Factorization mode                                             Default: PluFactLU              IN  */
    IPARM_STATIC_PIVOTING,       /**< Static pivoting                                                Default: -                         OUT */
    IPARM_FREE_CSCUSER,          /**< Free user CSC                                                  Default: 0                         IN  */
    IPARM_SCHUR_FACT_MODE,       /**< Specify if the Schur is factorized (@see plu_fact_mode_t)   Default: PluFactModeLocal       IN  */

    /* Solve */
    IPARM_SCHUR_SOLV_MODE,       /**< Specify the solve parts to apply (@see plu_solv_mode_t)     Default: PluSolvModeLocal       IN  */
    IPARM_APPLYPERM_WS,          /**< Enable/disable extra workspace for a thread-safe swap          Default: 1                         IN  */

    /* Refinement */
    IPARM_REFINEMENT,            /**< Refinement mode                                                Default: PluRefineGMRES         IN  */
    IPARM_NBITER,                /**< Number of iterations performed in refinement                   Default: -                         OUT */
    IPARM_ITERMAX,               /**< Maximum iteration number for refinement                        Default: 250                       IN  */
    IPARM_GMRES_IM,              /**< GMRES restart parameter                                        Default: 25                        IN  */

    /* Context */
    IPARM_SCHEDULER,             /**< Scheduler mode                                                 Default: PluSchedStatic         IN  */
    IPARM_THREAD_NBR,            /**< Number of threads per process (-1 for auto detect)             Default: -1                        IN  */
    IPARM_AUTOSPLIT_COMM,        /**< Automaticaly split communicator to have one MPI task by node   Default: 0                         IN  */

    /* Compression */
    IPARM_COMPRESS_MIN_WIDTH,    /**< Minimum width to compress a supernode                          Default: 120                       IN  */
    IPARM_COMPRESS_MIN_HEIGHT,   /**< Minimum height to compress an off-diagonal block               Default: 20                        IN  */
    IPARM_COMPRESS_WHEN,         /**< When to compress a supernode                                   Default: PluCompressNever       IN  */
    IPARM_COMPRESS_METHOD,       /**< Compression method (See plu_compress_method_t)              Default: PluCompressMethodPQRCP IN  */
    IPARM_COMPRESS_ORTHO,        /**< Orthogonalization method                                       Default: PluCompressOrthoCGS    IN  */
    IPARM_COMPRESS_RELTOL,       /**< Enable/Disable relative tolerance                              Default: 0                         IN  */

    /* MPI modes */
    IPARM_THREAD_COMM_MODE,      /**< Threaded communication mode                                    Default: PluThreadMultiple      IN  */

    /* Subset for old interface */
    IPARM_MODIFY_PARAMETER,      /**< Indicate if parameters have been set by user                   Default: 1                         IN  */
    IPARM_START_TASK,            /**< Indicate the first step to execute                             Default: PluTaskOrdering        IN  */
    IPARM_END_TASK,              /**< Indicate the last step to execute                              Default: PluTaskClean           IN  */
    IPARM_FLOAT,                 /**< Indicate the arithmetics                                       Default: PluDouble              IN  */
    IPARM_MTX_TYPE,              /**< Indicate matrix format                                         Default: -1                        IN  */
    IPARM_DOF_NBR,               /**< Degree of freedom per node                                     Default: 1                         IN  */
    IPARM_SIZE
} plu_iparm_t;

/**
 * @brief Float parameters
 */
typedef enum plu_dparm_e {
    DPARM_FILL_IN,               /**< Maximum memory (-DMEMORY_USAGE)                   Default: -                OUT */
    DPARM_EPSILON_REFINEMENT,    /**< Epsilon for refinement                            Default: -1.              IN  */
    DPARM_RELATIVE_ERROR,        /**< Relative backward error                           Default: -                OUT */
    DPARM_EPSILON_MAGN_CTRL,     /**< Epsilon for magnitude control                     Default: 0.               IN  */
    DPARM_ANALYZE_TIME,          /**< Time for Analyse step (wallclock)                 Default: -                OUT */
    DPARM_PRED_FACT_TIME,        /**< Predicted factorization time                      Default: -                OUT */
    DPARM_FACT_TIME,             /**< Time for Numerical Factorization step (wallclock) Default: -                OUT */
    DPARM_SOLV_TIME,             /**< Time for Solve step (wallclock)                   Default: -                OUT */
    DPARM_FACT_FLOPS,            /**< Factorization GFlops/s                            Default: -                OUT */
    DPARM_FACT_THFLOPS,          /**< Factorization theoretical Flops                   Default: -                OUT */
    DPARM_FACT_RLFLOPS,          /**< Factorization performed Flops                     Default: -                OUT */
    DPARM_SOLV_FLOPS,            /**< Solve GFlops/s                                    Default: -                OUT */
    DPARM_SOLV_THFLOPS,          /**< Solve theoretical Flops                           Default: -                OUT */
    DPARM_SOLV_RLFLOPS,          /**< Solve performed Flops                             Default: -                OUT */
    DPARM_REFINE_TIME,           /**< Time for Refinement step (wallclock)              Default: -                OUT */
    DPARM_A_NORM,                /**< ||A||_f norm                                      Default: -                OUT */
    DPARM_COMPRESS_TOLERANCE,    /**< Tolerance for low-rank kernels                    Default: 0.01             IN  */
    DPARM_COMPRESS_MIN_RATIO,    /**< Min ratio for rank w.r.t. strict rank             Default: 1.0              IN  */
    DPARM_SIZE
} plu_dparm_t;

/**
 * @brief Main steps for the plu() interface.
 *
 * Those enums are used of the IPARM_START_TASK and IPARM_END_TASK parameters
 * that configure the plu() call.
 */
typedef enum plu_task_e {
    PluTaskInit       = 0, /**< Startup the library          */
    PluTaskOrdering   = 1, /**< Ordering                     */
    PluTaskSymbfact   = 2, /**< Symbolic factorization       */
    PluTaskAnalyze    = 3, /**< Tasks mapping and scheduling */
    PluTaskNumfact    = 4, /**< Numerical factorization      */
    PluTaskSolve      = 5, /**< Numerical solve              */
    PluTaskRefine     = 6, /**< Numerical refinement         */
    PluTaskClean      = 7  /**< Clean                        */
} plu_task_t;

/**
 * @brief Verbose modes
 */
typedef enum plu_verbose_e {
    PluVerboseNot = 0, /**< Nothing  */
    PluVerboseNo  = 1, /**< Default  */
    PluVerboseYes = 2  /**< Extended */
} plu_verbose_t;

/**
 * @brief IO strategy for graph and ordering
 */
typedef enum plu_io_e {
    PluIONo         = 0, /**< No output or input */
    PluIOLoad       = 1, /**< Load ordering and symbol matrix instead of applying symbolic factorization step */
    PluIOSave       = 2, /**< Save ordering and symbol matrix after symbolic factorization step */
    PluIOLoadGraph  = 4, /**< Load graph  during ordering step */
    PluIOSaveGraph  = 8, /**< Save graph  during ordering step */
    PluIOLoadCSC    = 16,/**< Load CSC(d) during ordering step */
    PluIOSaveCSC    = 32 /**< Save CSC(d) during ordering step */
} plu_io_t;

/**
 * @brief Factorization Schur modes
 *
 * Describe which part of the matrix is factorized or not
 *
 */
typedef enum plu_fact_mode_e {
    PluFactModeLocal   = 0,
    PluFactModeSchur   = 1,
    PluFactModeBoth    = 2
} plu_fact_mode_t;

/**
 * @brief Solve Schur modes
 *
 * Describe which part of the solve is applied with the matrix
 *
 * \f[ A = \left( \begin{array}{cc}
 *             L_{11}U_{11} & U_{12} \\
 *             L_{21}       & S_{22} \end{array} \right) \f]
 *
 * For the lower part (and symmetrically for upper part):
 *   -# Solve \f[ L_{11} * x_{11} = b_{11} \f]
 *   -# Apply the update \f[ b_{22} = b_{22} - L_{21} * b_{11} \f]
 *   -# Solve the lower part of \f[ S_{22} * x_{22} = b_{22} \f] if S22 has been previously factorized.
 *
 * PluSolvModeLocal applies only the step 1.
 * PluSolvModeInterface applies steps 1 and 2.
 * PluSolvModeSchur applies all steps.
 *
 */
typedef enum plu_solv_mode_e {
    PluSolvModeLocal     = 0,
    PluSolvModeInterface = 1,
    PluSolvModeSchur     = 2
} plu_solv_mode_t;

/**
 * @brief Iterative refinement algorithms
 */
typedef enum plu_refine_e {
    PluRefineGMRES,   /**< GMRES              */
    PluRefineCG,      /**< Conjugate Gradient */
    PluRefineSR,      /**< Simple refinement  */
    PluRefineBiCGSTAB /**< BiCGStab           */
} plu_refine_t;

/**
 * @brief Arithmetic types.
 * This describes the different arithmetics that can be stored in a sparse matrix.
 * @remark The values start at 2 for compatibility purpose with PLASMA and
 * DPLASMA libraries, and they match the ones used in spm.
 * @sa spm_coeftype_t
 */
#define plu_coeftype_t spm_coeftype_t
#define PluPattern   SpmPattern
#define PluDouble    SpmDouble

/**
 * @brief Factorization algorithms available for IPARM_FACTORIZATION parameter
 */
typedef enum plu_factotype_e {
    PluFactGETRF = 2, /**< LU factorization                         */
    PluFactLU   = 2, /**< LU factorization                         */
} plu_factotype_t;

/**
 * @brief Scheduler
 */
typedef enum plu_scheduler_e {
    PluSchedSequential = 0, /**< Sequential                           */
    PluSchedStatic     = 1, /**< Shared memory with static scheduler  */
    PluSchedDynamic    = 4, /**< Shared memory with dynamic scheduler */
} plu_scheduler_t;

/**
 * @brief Ordering strategy
 */
enum plu_order_e {
    PluOrderMetis,    /**< Use Metis ordering                          */
    PluOrderPersonal, /**< Apply user's permutation, or load from file */
};

/**
 * @brief Error codes
 */
typedef enum plu_error_e {
    PLU_SUCCESS            = 0,  /**< No error                     */
    PLU_ERR_UNKNOWN        = 1,  /**< Unknown error                */
    PLU_ERR_ALLOC          = 2,  /**< Allocation error             */
    PLU_ERR_NOTIMPLEMENTED = 3,  /**< Not implemented feature      */
    PLU_ERR_OUTOFMEMORY    = 4,  /**< Not enough memory            */
    PLU_ERR_THREAD         = 5,  /**< Error with threads           */
    PLU_ERR_INTERNAL       = 6,  /**< Internal error               */
    PLU_ERR_BADPARAMETER   = 7,  /**< Bad parameters given         */
    PLU_ERR_FILE           = 8,  /**< Error in In/Out operations   */
    PLU_ERR_INTEGER_TYPE   = 9,  /**< Error with integer types     */
    PLU_ERR_IO             = 10, /**< Error with input/output      */
    PLU_ERR_MPI            = 11  /**< Error with MPI calls         */
} plu_error_t;

/**
 * @brief Compression strategy available for IPARM_COMPRESS_WHEN parameter
 */
typedef enum plu_compress_when_e {
    PluCompressNever,      /**< Do not use compression                                              */
    PluCompressWhenBegin,  /**< Compress before any numerical operation (Minimal-Memory)            */
    PluCompressWhenEnd,    /**< Compress after contributions were accumulated (Just-In-Time)        */
    PluCompressWhenDuring  /**< Compress after contributions from other supernodes were accumulated */
} plu_compress_when_t;

/**
 * @brief Compression method available for IPARM_COMPRESS_METHOD parameter
 */
typedef enum plu_compress_method_e {
    PluCompressMethodSVD,   /**< Use singular value decomposition for low-rank compression       */
    PluCompressMethodPQRCP, /**< Use partial QR with column pivoting for low-rank compression    */
    PluCompressMethodRQRCP, /**< Use randomized QR with column pivoting for low-rank compression */
    PluCompressMethodTQRCP, /**< Use truncated QR with column pivotingfor low-rank compression   */
    PluCompressMethodRQRRT, /**< Use randomized QR with rotation for low-rank compression        */
    PluCompressMethodNbr    /**< Total number of available compression methods                   */
} plu_compress_method_t;

/**
 * @brief Orthogonalization method available for IPARM_COMPRESS_ORTHO parameter
 */
typedef enum plu_compress_ortho_e {
    PluCompressOrthoCGS,        /**< Orthogonalize low-rank bases with Gram-Schimdt                                           */
    PluCompressOrthoQR,         /**< Orthogonalize low-rank bases with QR decomposition                                       */
    PluCompressOrthoPartialQR,  /**< Orthogonalize low-rank bases with projections in orthogonal space followed by smaller QR */
} plu_compress_ortho_t;

/**
 * @brief Splitting strategy available for IPARM_SPLITTING_STRATEGY parameter
 */
typedef enum plu_split_e {
    PluSplitNot,              /**< Do not apply dedicated low-rank clustering strategy */
    PluSplitKway,             /**< Use k-way partitioning                              */
    PluSplitKwayProjections,  /**< Use projections and k-way in clusters               */
} plu_split_t;

/**
 * @name Constants compatible with CBLAS & LAPACK & PLASMA
 *    The naming and numbering of the following constants is consistent with:
 *       - CBLAS from Netlib (http://www.netlib.org/blas/blast-forum/cblas.tgz)
 *       - C Interface to LAPACK from Netlib (http://www.netlib.org/lapack/lapwrapc/)
 *       - Plasma (http://icl.cs.utk.edu/plasma/index.html)
 */

/**
 * @brief Direction of the matrix storage
 */
typedef enum plu_layout_e {
    PluRowMajor  = 101, /**< Storage in row major order    */
    PluColMajor  = 102  /**< Storage in column major order */
} plu_layout_t;

/**
 * @brief Transpostion
 */
typedef enum plu_trans_e {
    PluNoTrans   = 111, /**< Use A         */
    PluTrans     = 112, /**< Use A^t       */
    PluConjTrans = 113  /**< Use conj(A^t) */
} plu_trans_t;

/**
 * @brief Matrix symmetry type property.
 * @remark Must match transposition.
 */
typedef enum plu_mtxtype_e {
    PluGeneral   = PluNoTrans,    /**< The matrix is general   */
    PluSymmetric = PluTrans,      /**< The matrix is symmetric */
    PluHermitian = PluConjTrans   /**< The matrix is hermitian */
} plu_mtxtype_t;

/**
 * @brief Upper/Lower part
 */
typedef enum plu_uplo_e {
    PluUpper      = 121, /**< Use lower triangle of A */
    PluLower      = 122, /**< Use upper triangle of A */
    PluUpperLower = 123  /**< Use the full A          */
} plu_uplo_t;

/**
 * @brief Data blocks used in the kernel
 */
typedef enum plu_coefside_e {
    PluLCoef  = 0, /**< Coefficients of the lower triangular L are used         */
    PluUCoef  = 1, /**< Coefficients of the upper triangular U are used         */
    PluLUCoef = 2  /**< Coefficients of the upper/lower triangular U/L are used */
} plu_coefside_t;

/**
 * @brief Diagonal
 */
typedef enum plu_diag_e {
    PluNonUnit = 131, /**< Diagonal is non unitary */
    PluUnit    = 132  /**< Diagonal is unitary     */
} plu_diag_t;

/**
 * @brief Side of the operation
 */
typedef enum plu_side_e {
    PluLeft  = 141, /**< Apply operator on the left  */
    PluRight = 142  /**< Apply operator on the right */
} plu_side_t;

/**
 * @brief Norms
 */
typedef enum plu_normtype_e {
    PluOneNorm       = 171, /**< One norm:       max_j( sum_i( |a_{ij}| ) )   */
    PluFrobeniusNorm = 174, /**< Frobenius norm: sqrt( sum_{i,j} (a_{ij}^2) ) */
    PluInfNorm       = 175, /**< Inifinite norm: max_i( sum_j( |a_{ij}| ) )   */
    PluMaxNorm       = 177  /**< Inifinite norm: max_{i,j}( | a_{ij} | )      */
} plu_normtype_t;

/**
 * @brief Direction
 */
typedef enum plu_dir_e {
    PluDirForward  = 391, /**< Forward direction   */
    PluDirBackward = 392, /**< Backward direction  */
} plu_dir_t;

#endif /* _plu_api_h_ */