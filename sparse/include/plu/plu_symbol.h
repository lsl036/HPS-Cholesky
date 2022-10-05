
#ifndef _isched_barrier_h_
#define _isched_barrier_h_

#include <unistd.h>
#include <pthread.h>

/* The Linux includes are completely screwed up right now. Even if they
 * correctly export a _POSIX_BARRIER define the barrier functions are
 * not correctly defined in the pthread.h. So until we figure out
 * how to correctly identify their availability, we will have to
 * disable them.
 */
#if defined(_POSIX_BARRIERS) && (_POSIX_BARRIERS - 20012L) >= 0 && 0

BEGIN_C_DECLS

typedef pthread_barrier_t isched_barrier_t;
#define isched_barrier_init      pthread_barrier_init
#define isched_barrier_wait      pthread_barrier_wait
#define isched_barrier_destroy   pthread_barrier_destroy
#define ISCHED_IMPLEMENT_BARRIERS 0

#else

typedef struct isched_barrier_s {
    int                 count;
    volatile int        curcount;
    volatile int        generation;
    pthread_mutex_t     mutex;
    pthread_cond_t      cond;
} isched_barrier_t;

int isched_barrier_init(isched_barrier_t *barrier, const void *pthread_mutex_attr, unsigned int count);
int isched_barrier_wait(isched_barrier_t*);
int isched_barrier_destroy(isched_barrier_t*);
#define ISCHED_IMPLEMENT_BARRIERS 1

END_C_DECLS

#endif


#endif /* _isched_barrier_h_ */


#ifndef _isched_h_
#define _isched_h_

BEGIN_C_DECLS

enum isched_action_e {
    ISCHED_ACT_STAND_BY,
    ISCHED_ACT_PARALLEL,
    ISCHED_ACT_FINALIZE
};

struct isched_s;
typedef struct isched_s isched_t;

/**
 * Thread structure of the execution context of one instance of the scheduler
 */
typedef struct isched_thread_s {
    isched_t        *global_ctx;
    int              rank;
    int              bindto;
} isched_thread_t;

/**
 * Global structure of the execution context of one instance of the scheduler
 */
struct isched_s {
    int              world_size;

    isched_barrier_t barrier;
    pthread_mutex_t  statuslock;
    pthread_cond_t   statuscond;
    volatile int     status;

    pthread_t       *tids;
    isched_thread_t *master;

    void           (*pfunc)(isched_thread_t*, void*);
    void            *pargs;
};

#define isched_topo_init               isched_nohwloc_init
#define isched_topo_destroy            isched_nohwloc_destroy
#define isched_topo_bind_on_core_index isched_nohwloc_bind_on_core_index
#define isched_topo_unbind             isched_nohwloc_unbind
#define isched_topo_world_size         isched_nohwloc_world_size

int  isched_topo_init(void);
int  isched_topo_destroy(void);
int  isched_topo_bind_on_core_index(int);   
int  isched_topo_unbind();
int  isched_topo_world_size();

static inline void
isched_parallel_call( isched_t *isched, void (*func)(isched_thread_t*, void*), void *args )
{
    pthread_mutex_lock(&isched->statuslock);
    isched->pfunc  = func;
    isched->pargs  = args;
    isched->status = ISCHED_ACT_PARALLEL;
    pthread_mutex_unlock(&isched->statuslock);
    /* Wake up all threads waiting for condition variables COND.  */
    pthread_cond_broadcast(&isched->statuscond);
    isched_barrier_wait( &(isched->barrier) );
    isched->status = ISCHED_ACT_STAND_BY;
    func( isched->master, args );
    isched_barrier_wait( &(isched->barrier) );
}

isched_t *ischedInit(int cores, const int *coresbind);
int ischedFinalize(isched_t *isched);

END_C_DECLS

#endif /* _isched_h_ */


#ifndef _symbol_h_
#define _symbol_h_

/**
 * @brief Symbol column block structure.
 */
typedef struct symbol_cblk_s {
    Sparse_long fcolnum; /**< First column index               */
    Sparse_long lcolnum; /**< Last column index (inclusive)    */
    Sparse_long bloknum; /**< First block in column (diagonal) */
    Sparse_long brownum; /**< First block in row facing the diagonal block in browtab, 0-based */
    int8_t       selevtx;
#if defined( PLU_SYMBOL_DUMP_SYMBMTX )
    Sparse_long split_cblk;
#endif
} symbol_cblk_t;

/**
 * @brief Symbol block structure.
 */
typedef struct symbol_blok_s {
    Sparse_long frownum; /**< First row index            */
    Sparse_long lrownum; /**< Last row index (inclusive) */
    Sparse_long lcblknm; /**< Local column block         */
    Sparse_long fcblknm; /**< Facing column block        */
} symbol_blok_t;

/**
 * @brief Symbol matrix structure.
 *
 * This structure describes the symbolic block structure of the factorized
 * matrix L, U is never stored as it is a symmetry of L. This structure is
 * global and replicated on all processes. The default way to number the block
 * is the CSC format where block are continuously number per column, the browtab
 * array stores the CSR representation of the L structure to have a faster
 * access to the list of blocks updating a column block.
 *
 */
typedef struct symbol_matrix_s {
    Sparse_long   baseval;   /**< Base value for numbering                   */
    Sparse_long   cblknbr;   /**< Number of column blocks                    */
    Sparse_long   bloknbr;   /**< Number of blocks                           */
    Sparse_long   nodenbr;   /**< Number of nodes (Equal to gN in spm)       */
    Sparse_long   schurfcol; /**< First column of the schur complement       */
    symbol_cblk_t *cblktab;   /**< Array of column blocks [+1,based]          */
    symbol_blok_t *bloktab;   /**< Array of blocks in CSC format [based]      */
    Sparse_long  *browtab;   /**< Array of blocks in CSR format [based]      */
    Sparse_long   dof;       /**< Degrees of freedom per node (constant
                                   if > 0, variadic if < 1                    */
    Sparse_long  *dofs;      /**< Array of the first column of each element
                                   in the expanded matrix [+1,based]          */
} symbol_matrix_t;

/**
 * @name Symbol basic subroutines
 * @{
 */
void pluSymbolInit   ( const plu_graph_t  *graph,
                          const plu_order_t  *order,
                                symbol_matrix_t *symbptr );
void pluSymbolExit   (       symbol_matrix_t *symbptr );
void pluSymbolBase   (       symbol_matrix_t *symbptr,
                          const Sparse_long     baseval );
void pluSymbolRealloc(       symbol_matrix_t *symbptr );
int  pluSymbolCheck  ( const symbol_matrix_t *symbptr );
void pluSymbolExpand (       symbol_matrix_t *symbptr );

/**
 * @}
 * @name Symbol IO subroutines
 * @{
 */
int pluSymbolSave( const symbol_matrix_t *symbptr, FILE *stream );
int pluSymbolLoad(       symbol_matrix_t *symbptr, FILE *stream );
int pluSymbolDraw( const symbol_matrix_t *symbptr, FILE *stream );

/**
 * @}
 * @name Symbol statistical information subroutines
 * @{
 */
void         pluSymbolPrintStats( const symbol_matrix_t *symbptr );
Sparse_long pluSymbolGetNNZ( const symbol_matrix_t *symbptr );
void         pluSymbolGetFlops( const symbol_matrix_t *symbmtx,
                                   plu_coeftype_t      flttype,
                                   plu_factotype_t     factotype,
                                   double                *thflops,
                                   double                *rlflops );
void         pluSymbolGetTimes( const symbol_matrix_t *symbmtx,
                                   plu_coeftype_t      flttype,
                                   plu_factotype_t     factotype,
                                   double                *cblkcost,
                                   double                *blokcost );

/**
 * @}
 * @name Symbol reordering subroutines
 * @{
 */
void pluSymbolReordering( plu_data_t * );
void pluSymbolReorderingPrintComplexity( const symbol_matrix_t *symbptr );

/**
 * @}
 * @name Symbol construction subroutines
 * @{
 */
int          pluSymbolFaxDirect ( symbol_matrix_t      *symbptr,
                                     const plu_graph_t *graphA,
                                     const plu_order_t *ordeptr );
int          pluSymbolFaxILUk   ( symbol_matrix_t      *symbptr,
                                     Sparse_long          levelk,
                                     const plu_graph_t *graphA,
                                     const plu_order_t *ordeptr );
void         pluSymbolRustine   ( symbol_matrix_t *symbptr, symbol_matrix_t *symbptr2 );
void         pluSymbolBuildRowtab( symbol_matrix_t *symbptr );
Sparse_long pluSymbolGetFacingBloknum( const symbol_matrix_t *symbptr,
                                           Sparse_long           bloksrc,
                                           Sparse_long           bloknum,
                                           Sparse_long           startsearch,
                                           int                    ricar );

#endif /* _symbol_h_ */

#ifndef _queue_h_
#define _queue_h_

/**
 * @brief Queue item structure.
 */
typedef struct plu_queue_item_s {
    double       key1;   /**< Key 1 of the element   */
    double       key2;   /**< Key 2 of the element   */
    Sparse_long eltptr; /**< Pointer to the element */
} plu_queue_item_t;

/**
 * @brief Queue structure.
 */
typedef struct plu_queue_s {
    Sparse_long          size;   /**< Allocated memory size                           */
    volatile Sparse_long used;   /**< Number of element in the queue                  */
    plu_queue_item_t  *elttab; /**< Array of the element                            */
    plu_atomic_lock_t  lock;   /**< Lock for insertion and removal in shared memory */
} plu_queue_t;

int          pqueueInit(        plu_queue_t *, Sparse_long );
void         pqueueExit(        plu_queue_t * );
Sparse_long pqueueSize(  const plu_queue_t * );
void         pqueueClear(       plu_queue_t * );
void         pqueuePush2(       plu_queue_t *, Sparse_long, double, double );
Sparse_long pqueueRead ( const plu_queue_t * );
Sparse_long pqueuePop2 (       plu_queue_t *, double *, double * );
void         pqueuePrint( const plu_queue_t * );

/**
 * @brief Push an element with a single key.
 * @param[inout] q
 *               The queue structure.
 * @param[in]    elt
 *               The element to insert.
 * @param[in]    key1
 *               The first key of the element to insert (the second will be 0.).
 */
static inline void
pqueuePush1(plu_queue_t *q, Sparse_long elt, double key1) {
    pqueuePush2( q, elt, key1, 0. );
}

/**
 * @brief Pop the head of the queue whithout returning the keys.
 * @param[inout] q
 *               The queue structure.
 * @return The element at the head of the queue.
 */
static inline Sparse_long
pqueuePop(plu_queue_t *q){
    return pqueuePop2(q, NULL, NULL);
}

/**
 * @brief Pop the head of the queue and get the associated first key.
 * @param[inout] q
 *               The queue structure.
 * @param[out]   key1
 *               The first key of the element removed from the head of the queue.
 * @return The element at the head of the queue.
 */
static inline Sparse_long
pqueuePop1(plu_queue_t *q, double *key1){
    return pqueuePop2(q, key1, NULL);
}

#endif /* _queue_h_ */

#ifndef _pludata_h_
#define _pludata_h_

/*
 * Steps of the plu solver
 */
#define STEP_INIT      (1 << 0)
#define STEP_ORDERING  (1 << 1)
#define STEP_SYMBFACT  (1 << 2)
#define STEP_ANALYSE   (1 << 3)
#define STEP_CSC2BCSC  (1 << 4)
#define STEP_BCSC2CTAB (1 << 5)
#define STEP_NUMFACT   (1 << 6)
#define STEP_SOLVE     (1 << 7)
#define STEP_REFINE    (1 << 8)

struct plu_bcsc_s;
typedef struct plu_bcsc_s plu_bcsc_t;

struct plu_model_s;
typedef struct plu_model_s plu_model_t;

/**
 * @brief Main PLU data structure
 * This structure holds all informations related to the library and problem
 * instance. It stores information from one step to another.
 * @warning This structure should not be modified directly by the user.
 *
 */
struct plu_data_s {
    Sparse_long    *iparm;              /**< Store integer parameters (input/output)                             */
    double          *dparm;              /**< Store floating parameters (input/output)                            */

    Sparse_long     steps;              /**< Bitmask of the steps performed or not                               */

    isched_t        *isched;             /**< Internal scheduler structure that is always available               */
    const sparse_csc *csc;               /**< Pointer to the user csc structure used as input                     */

    plu_graph_t  *graph;              /**< Symmetrized graph of the problem used within ordering
                                              and symbolic factorization steps.                                   */
    Sparse_long     schur_n;            /**< Number of entries for the Schur complement                          */
    Sparse_long    *schur_list;         /**< List of entries for the schur complement                            */
    Sparse_long     zeros_n;            /**< Number of diagonal entries considered as zeros                      */
    Sparse_long    *zeros_list;         /**< List of diagonal entries considered as zeros                        */
    plu_order_t  *ordemesh;           /**< Ordering structure                                                  */

    symbol_matrix_t *symbmtx;            /**< Symbol Matrix                                                       */

    plu_bcsc_t   *bcsc;               /**< Csc after reordering grouped by cblk                                */
    SolverMatrix    *solvmatr;           /**< Solver informations associted to the matrix problem                 */

    plu_model_t  *cpu_models;         /**< CPU model coefficients for the kernels                              */
    
    char            *dir_global;         /**< Unique directory name to store output files                         */
    char            *dir_local;          /**< Unique directory name to store output specific to a MPI process     */

    /* Backup for old plu interface */
    void            *b;
    void            *x0;

    /**
     * Former fields that are no longer used for now
     */

    int             *bindtab;            /*+ Tabular giving for each thread a CPU to bind it too                 */
    void            *schur_tab;
    Sparse_long     schur_tab_set;
    int              cscInternFilled;
    int              scaling;            /*+ Indicates if the matrix has been scaled                             */
    void            *scalerowtab;        /*+ Describes how the matrix has been scaled                            */
    void            *iscalerowtab;
    void            *scalecoltab;
    void            *iscalecoltab;
#ifdef WITH_SEM_BARRIER
    sem_t           *sem_barrier;        /*+ Semaphore used for AUTOSPLIT_COMM barrier                           */
#endif
    Sparse_long     plu_id;          /*+ Id of the plu instance (PID of first MPI task)                   */
};

#endif /* _pludata_h_ */


#ifndef _symbol_cost_h_
#define _symbol_cost_h_

/**
 * @brief Cost functions to compute statistics on the symbolic structure
 */
typedef struct symbol_function_s {
    double (*diag     )(Sparse_long);                             /**< Return a statistic based on the diagonal block */
    double (*trsm     )(Sparse_long, Sparse_long);               /**< Return a statistic based on the sum of all
                                                                        off-diagonal of each column-block              */
    double (*update   )(Sparse_long, Sparse_long);               /**< Return a statistic for a large accumulated
                                                                        update per column-block                        */
    double (*blkupdate)(Sparse_long, Sparse_long, Sparse_long); /**< Return a statistic for each individual
                                                                        off-diagonal block                             */
} symbol_function_t;

/**
 * @brief array of pointer to the flops functions per factorization and arithmetic
 */
extern symbol_function_t flopstable;
/**
 * @brief array of pointer to the performance functions per factorization and arithmetic
 */
extern symbol_function_t perfstable;

#endif /* _symbol_cost_h_ */

#ifndef _symbol_fax_h_
#define _symbol_fax_h_

/**
 * @brief Prime number for hashing vertex numbers.
 */
#define SYMBOL_FAX_HASHPRIME        17

/**
 * @brief The chained column block structure.
 *
 * These blocks are chained in a single linked list
 * for block merge with blocks of left columns.
 *
 */
typedef struct symbol_faxtlok_s {
  Sparse_long frownum; /**< First row index            */
  Sparse_long lrownum; /**< Last row index (inclusive) */
  Sparse_long fcblknm; /**< Facing column block        */
  Sparse_long nextnum; /**< Index of next block        */
} SymbolFaxTlok;

/**
 * @}
 */
#endif /* _symbol_fax_h_ */

#ifndef _symbol_reorder_h_
#define _symbol_reorder_h_

void
symbol_reorder_cblk( const symbol_matrix_t *symbptr,
                     const symbol_cblk_t   *cblk,
                     plu_order_t        *order,
                     const Sparse_long    *levels,
                     Sparse_long          *depthweight,
                     Sparse_long           depthmax,
                     Sparse_long           split_level,
                     Sparse_long           stop_criterion );

void
symbol_reorder( plu_data_t *plu_data,
                Sparse_long   maxdepth,
                Sparse_long  *levels );

#endif /* _symbol_reorder_h_ */

#ifndef _elimintree_h_
#define _elimintree_h_

/**
 * @brief Node of the elimination tree.
 */
typedef struct etree_node_s {
    double       total;   /**< Cost of the treenode only (compute + send)    */
    double       subtree; /**< Cost of the subtree (includes total)          */
    double       cripath; /**< Cost of the citical path to the node included */
    int          ndlevel; /**< Node depth in the elimination tree            */
    int          sonsnbr; /**< Number of sons                                */
    Sparse_long fathnum; /**< index of the father node                      */
    Sparse_long fsonnum; /**< index of first son                            */
} eTreeNode_t;

/**
 * @brief Elimination tree.
 */
typedef struct etree_s {
    Sparse_long   baseval; /**< Base value for numberings         */
    Sparse_long   nodenbr; /**< Number of nodes                   */
    eTreeNode_t  * nodetab; /**< Array of node          [+1,based] */
    Sparse_long * sonstab; /**< Sons index of nodes               */
} EliminTree;

EliminTree   *eTreeInit      (      Sparse_long);
void          eTreeExit      (      EliminTree *);
void          eTreeGenDot    (const EliminTree *, FILE *);
void          eTreePrint     (const EliminTree *, FILE *, Sparse_long );
void          eTreeSetSons   (      EliminTree *);
Sparse_long  eTreeLeavesNbr (const EliminTree *);
Sparse_long  eTreeLevel     (const EliminTree *);
Sparse_long  eTreeNodeLevel (const EliminTree *, Sparse_long );
EliminTree   *eTreeBuild     (const symbol_matrix_t *);

Sparse_long eTreeComputeLevels   ( EliminTree *, Sparse_long, Sparse_long );
Sparse_long eTreeGetLevelMinIdx  ( const EliminTree *, Sparse_long, Sparse_long, Sparse_long );

/**
 * @brief Return the father of a given node.
 **/
static inline eTreeNode_t *
eTreeFather( const EliminTree *etree, Sparse_long node )
{
    return etree->nodetab + etree->nodetab[node].fathnum;
}

/**
 * @brief Return the i^{th} son of a given node.
 **/
static inline Sparse_long
eTreeSonI( const EliminTree *etree, Sparse_long node, Sparse_long i )
{
    return etree->sonstab[ etree->nodetab[node].fsonnum + i ];
}

/**
 * @brief Return the root of the elimination tree.
 **/
static inline Sparse_long
eTreeRoot( const EliminTree *etree )
{
    (void)etree;
    return -1;
}

#endif /* _elimintree_h_ */


#ifndef _extendvector_h_
#define _extendvector_h_

/**
 * @brief The extend integer array structure.
*/
typedef struct ExtendVectorINT_s {
    Sparse_long  vecsize; /**< The size of the vector             */
    Sparse_long  eltnbr;  /**< The number of elements stored      */
    Sparse_long *inttab;  /**< The actual array with the elements */
} ExtendVectorINT;

Sparse_long *extendint_Init  (       ExtendVectorINT *, Sparse_long );
void          extendint_Exit  (       ExtendVectorINT * );
void          extendint_Add   (       ExtendVectorINT *, Sparse_long );
Sparse_long  extendint_Size  ( const ExtendVectorINT * );
Sparse_long  extendint_Read  ( const ExtendVectorINT *, Sparse_long );
void          extendint_Clear (       ExtendVectorINT * );
void          extendint_ToSize(       ExtendVectorINT *, Sparse_long );
void          extendint_incr  (       ExtendVectorINT * );

#endif /* _extendvector_h_ */

#ifndef _simu_timer_h_
#define _simu_timer_h_

/**
 * @brief Timer for the simulation.
 */
typedef struct simu_timer_s {
    double s; /**< Second in the timer */
    /*  double ms;*/
} SimuTimer;

/**
 * @brief Compare two timings
 * @param[in] t1
 *            The first timer.
 * @param[in] t2
 *            The second timer.
 * @return True if the t1 is smaller than t2, false otherwise.
 */
static inline int
timerComp(const SimuTimer *t1,
          const SimuTimer *t2)
{
    /* Return (t1 < t2) */
    if(t1->s < t2->s) {
        return 1;
    }
    else {
        return 0;
    }
}

/**
 * @brief Increment the timer
 * @param[inout] timer
 *               The timer to update.
 * @param[in]    t
 *               The time to add to the timer.
 */
static inline void
timerAdd(SimuTimer *timer, double t)
{
    timer->s += t;
}

/**
 * @brief Get the timer value
 * @param[in] timer
 *            The timer to read.
 * @return The timer value in second.
 */
static inline double
timerVal(const SimuTimer *timer)
{
    return timer->s;
}

/**
 * @brief Set the timer value
 * @param[inout] timer
 *               The timer to set
 * @param[in]    t
 *               The value to set
 */
static inline void
timerSet(SimuTimer *timer, double t)
{
    timer->s = t;
}

/**
 * @brief Set the timer value if the value is greater than the actual one.
 * @param[inout] timer
 *               The timer to update
 * @param[in]    t
 *               The time to compare with and set if larger than the timer.
 */
static inline void
timerSetMax(SimuTimer *timer, double t)
{
    if ( t > timer->s ) {
        timer->s = t;
    }
}

#endif /* _simu_timer_h_ */

#ifndef _cost_h_
#define _cost_h_

/**
 * @brief Arrays of double to store the cost of each element in the matrix
 */
typedef struct cost_matrix_s {
    double *blokcost; /**< Cost of the update generated by this block
                           for off-diagonal block, fact+solve otherwise */
    double *cblkcost; /**< Cost of all the operations linked to a panel */
} CostMatrix;

void        costMatrixInit ( CostMatrix *costmtx );
void        costMatrixExit ( CostMatrix *costmtx );
CostMatrix *costMatrixBuild( const symbol_matrix_t *symbmtx,
                             plu_coeftype_t      flttype,
                             plu_factotype_t     factotype );

#endif /* _cost_h_ */

#ifndef _cand_h_
#define _cand_h_

/**
 * @brief Processor candidate group to own a column blok
 */
typedef struct cand_s {
    double       costlevel; /**< Cost from root to node                              */
    Sparse_long treelevel; /**< Level of the cblk in the elimination tree (depth from the root) */
    Sparse_long fcandnum;  /**< first processor number of this candidate group      */
    Sparse_long lcandnum;  /**< last processor number of this candidate group       */
    Sparse_long fccandnum; /**< first cluster number of the cluster candidate group */
    Sparse_long lccandnum; /**< last cluster number of the cluster candidate group  */
    Sparse_long cluster;   /**< Cluster id on which the task will be executed       */
    int8_t       cblktype;  /**< type of the distribution                            */
} Cand;

Cand *candInit (       Sparse_long     cblknbr );
void  candExit (       Cand            *candtab );
int   candCheck( const Cand            *candtab,
                 const symbol_matrix_t *symbmtx );
void  candSave ( const Cand            *candtab,
                       Sparse_long     cblknbr,
                 const char            *directory );
void  candBuild( Sparse_long           level_tasks2d,
                 Sparse_long           width_tasks2d,
                 plu_compress_when_t lr_when,
                 Sparse_long           lr_width,
                 Cand                  *candtab,
                 EliminTree            *etree,
                 const symbol_matrix_t *symbmtx,
                 const CostMatrix      *costmtx );


void candUpdate         ( Cand                  *candtab,
                          EliminTree            *etree,
                          const symbol_matrix_t *symbmtx,
                          const CostMatrix      *costmtx );

void candSetClusterCand(       Cand          *candtab,
                               Sparse_long   cblknbr,
                         const Sparse_long  *core2clust,
                               Sparse_long   coresnbr );

void candGenDot          ( const EliminTree *etree,
                           const Cand       *candtab,
                           FILE             *stream );
void candGenDotLevel     ( const EliminTree *etree,
                           const Cand       *candtab,
                           FILE             *stream,
                           Sparse_long      level );
void candGenCompressedDot( const EliminTree *etree,
                           const Cand       *candtab,
                           FILE             *stream );

#endif /* _cand_h_ */

#include "plu_solver.h"

#ifndef _simu_h_
#define _simu_h_

/**
 * @brief Process structure for the simulation.
 */
typedef struct simu_cluster_s {
    Sparse_long     fprocnum;   /**< Global index of the first processor belonging to the cluster (Check is it is not redundant) */
    Sparse_long     lprocnum;   /**< Global index of the last processor belonging to the cluster (Check is it is not redundant)  */
    ExtendVectorINT *ftgtsend;   /**< Arrays of ftgt sent by this proc (one vector per processor)                                 */
    Sparse_long     prionum;    /**< Counter to order tasks on one cluster                                                       */
} SimuCluster;

/**
 * @brief Thread structure for the simulation.
 */
typedef struct simu_proc_s {
    SimuTimer        timer;      /**< Simulated clock of the processor                                  */
    plu_queue_t  *readytask;  /**< Heap of tasks ready to be executed                                */
    plu_queue_t  *futuretask; /**< Heap of tasks ready to be executed in a near future (after timer) */
    ExtendVectorINT *tasktab;    /**< Vector to store tasks affected to the candidate                   */
    char            *procalias;  /**< Paje trace alias to the processor if PLU_BLEND_GENTRACE is eenabled */
} SimuProc;

/**
 * @brief Fan-in structure for the simulation.
 */
typedef struct simu_ftgt_s {
    solver_ftgt_t ftgt;         /**< Fan-in informations                            */
    Sparse_long  clustnum;     /**< Cluster sending the contribution               */
    SimuTimer     timerecv;     /**< Simulated clock of the reception time          */
    double        costsend;     /**< Cost to send the contribution                  */
    double        costadd;      /**< Cost to add the contribution to its final cblk */
} SimuFtgt;

/**
 * @brief Column block structure for the simulation.
 */
typedef struct simu_cblk_s {
    Sparse_long ctrbcnt;       /**< Counter of remaining contributions for the cblk */
} SimuCblk;

/**
 * @brief Block structure for the simulation.
 */
typedef struct simu_blok_s {
    Sparse_long tasknum;       /**< Task index opeating on this block (stored per block for 2D computations)   */
    Sparse_long ftgtnum;       /**< Index of the first fanin destinated to this
                                     block in the ftgttab. This index is also used to find the first cblk timer
                                     (one per cand proc) in the timetab array                                   */
    Sparse_long ctrbcnt;       /**< Counter of remaining contributions                                         */
    int          fccandnum;     /**< First candidate that is attributed to the cblk of the block                */
    int          ownerclust;    /**< Processor on which the block is distributed                                */
} SimuBlok;

/**
 * @brief Task structure for the simulation.
 */
typedef struct simu_task_s {
    Sparse_long prionum;       /**< priority of the task                                      */
    Sparse_long cblknum;       /**< Number of the cblknum the task deals with                 */
    Sparse_long bloknum;       /**< number of the block that the task deals with              */
    Sparse_long bloknum2;      /**< */
    Sparse_long facebloknum;   /**< Number of the facing block for E2                         */
    SimuTimer    time;          /**< Time the task is ready if it doesn't need message         */
    Sparse_long mesglen;       /**< Time to send the block target                             */
    double       cost;          /**< Cost of the task                                          */
    Sparse_long ctrbcnt;       /**< nbr ftgt + le btgt (+ E1 pret si E2)                      */
    Sparse_long ftgtcnt;       /**< nbr of contrib from fan-in target                         */
    Sparse_long tasknext;      /**< chainage des E1 ou E2, si fin = -1 => liberer les btagptr */
} SimuTask;

/**
 * @brief Control structure for the simulation.
 */
typedef struct simuctrl_s {
    Sparse_long  cblknbr;      /**< Number of cblk                                            */
    Sparse_long  ftgtprio;     /**< Priority to assign to current ftgts                       */
    Sparse_long  tasknbr;      /**< Number of tasks                                           */
    Sparse_long  ftgtcnt;      /**< Number of received communication                          */
    SimuTask     *tasktab;      /**< SimuTask vector                                           */
    SimuProc     *proctab;      /**< Virtual processor tab                                     */
    SimuCluster  *clustab;      /**< Virtual cluster tab                                       */
    Sparse_long *ownetab;      /**< Vector containing the distribution of the diagonal blok   */
    SimuCblk     *cblktab;      /**< SimuCblk vector                                           */
    SimuBlok     *bloktab;      /**< SimuBlok vector                                           */
    SimuFtgt     *ftgttab;      /**< Vector containing the fan in target                       */
    Sparse_long  ftgtnbr;      /**< The number of fan-in contribution                         */
    SimuTimer    *ftgttimetab;  /**< Vector containing a timer for each cluster on each ftgt   */
} SimuCtrl;


Sparse_long simuInit        ( SimuCtrl *, const symbol_matrix_t *, const Cand *, Sparse_long, Sparse_long );
Sparse_long simuRealloc     ( SimuCtrl *, Sparse_long, Sparse_long );
void         simuExit        ( SimuCtrl *, Sparse_long, Sparse_long, Sparse_long );
void         simuTaskBuild   ( SimuCtrl *, const symbol_matrix_t * );

#define CLUST2INDEX(n,c) ((c) + simuctrl->bloktab[n].ftgtnum - simuctrl->bloktab[n].fccandnum)
#define INDEX2CLUST(r,s) ((r) - simuctrl->bloktab[s].ftgtnum + simuctrl->bloktab[s].fccandnum)
#define TIMER(pr)        (&(simuctrl->proctab[pr].timer))

#endif /* _simu_h_ */



#ifndef _blendctrl_h_
#define _blendctrl_h_

/**
 * @brief The type and structure definitions.
*/
typedef struct blendctrl_s {
    Sparse_long    count_ops ;      /**< Print costs in term of number of elementary operations            */
    Sparse_long    debug ;          /**< Make additional checks after each step                            */
    Sparse_long    timer;           /**< Print execution times                                             */
    Sparse_long    ooc;             /**< Enable the out-of-core version of Plu (Option unused for now)  */
    Sparse_long    ricar;           /**< Enable the ILU(k) dedicated steps                                 */
    Sparse_long    leader;          /**< Leader for sequential tasks                                       */

    Sparse_long    allcand;         /**< All processors are candidate for each cblk                        */
    Sparse_long    nocrossproc;     /**< Forbid a processor to be candidate in two
                                          different branches shared with different partners                 */
    Sparse_long    costlevel;       /**< Enable/disable computation and use of subtree cost                */

    Sparse_long    blcolmin ;       /**< Minimun number of columns for a good use of BLAS primitives       */
    Sparse_long    blcolmax;        /**< Maximum number of columns for a good use of BLAS primitives       */
    Sparse_long    abs;             /**< Adaptative block size:
                                            - 0, all block are cut to blcolmin
                                            - >0, try to make (ncand*abs) cblk                              */
    Sparse_long    up_after_split;  /**< Update the costmtx and candtab arrays after splitting the symbol  */

    Sparse_long    level_tasks2d;   /**< Level to shift from 1D to 2D. Automaticaly computed if < 0, only 1D if 0 */
    Sparse_long    width_tasks2d;   /**< Minimal width to consider a cblk 2D if autolevel (level_tasks2d < 0)     */

    Sparse_long    clustnum;        /**< Id of current MPI process                                         */
    Sparse_long    clustnbr;        /**< Number of MPI processes                                           */
    Sparse_long    total_nbcores;   /**< Total number of physical cores used for the simulation            */
    Sparse_long    total_nbthrds;   /**< Total number of threads used for the simulation                   */
    Sparse_long    local_nbcores;   /**< Local number of physical cores used by the current MPI process    */
    Sparse_long    local_nbthrds;   /**< Local number of threads used by the current MPI process           */
    Sparse_long    local_nbctxts;   /**< Local number of contexts (used for dynamic scheduler and runtimes)*/
    Sparse_long   *clust2smp;       /**< clust2smp[i] = SMP node on which i_th MPI
                                          process is running, if multiple MPI processes per node            */
    Sparse_long   *core2clust;      /**< core2clust[i] = cluster owning the core i                         */

    Sparse_long   *iparm;           /**< In/Out Integer parameters                                         */
    double         *dparm;           /**< In/Out Float parameters                                           */
    const char     *dirname;         /**< Temporary unique directory to store output files
                                          (Initialized to plu_data->dir_local)                           */

    EliminTree        *etree;        /**< the elimination tree                                              */
    CostMatrix        *costmtx;      /**< the cost bounded to each cblk and blok                            */
    Cand              *candtab;      /**< processor candidate tab                                           */
    FILE              *tracefile;    /**< File holding the simulated trace                                  */

} BlendCtrl;

int  blendCtrlInit ( plu_data_t *plu_data,
                     BlendCtrl     *ctrl );

void blendCtrlExit (BlendCtrl *);

void getCommunicationCosts( const BlendCtrl *ctrl,
                            Sparse_long clustsrc,
                            Sparse_long clustdst,
                            Sparse_long sync_comm_nbr,
                            double *startup,
                            double *bandwidth);

#endif /* _blendctrl_h_ */

#ifndef _blend_h_
#define _blend_h_

void propMappTree   ( Cand             *candtab,
                      const EliminTree *etree,
                      Sparse_long      candnbr,
                      int               nocrossproc,
                      int               allcand );
void splitSymbol    ( BlendCtrl    *ctrl,
                      symbol_matrix_t *symbmtx );
void simuRun        ( SimuCtrl *,
                      const BlendCtrl *,
                      const symbol_matrix_t * );
#endif /* _blend_h_ */


#ifndef _fax_csr_h_
#define _fax_csr_h_

/**
 * @brief Fax blocked csr structure
 */
typedef struct fax_csr_s {
    Sparse_long   n;
    Sparse_long   total_nnz;
    Sparse_long * nnz;
    Sparse_long **rows;
} fax_csr_t;

void         faxCSRInit( Sparse_long n, fax_csr_t *csr );
void         faxCSRClean( fax_csr_t *csr );

Sparse_long faxCSRGetNNZ( const fax_csr_t *csr );

int  faxCSRGenPA( const plu_graph_t *graphA, const Sparse_long *perm, fax_csr_t *graphPA );
void faxCSRCompact( fax_csr_t *csr );

void faxCSRCblkCompress( const fax_csr_t      *graphA,
                         const plu_order_t *order,
                         fax_csr_t            *graphL,
                         Sparse_long         *work );

Sparse_long faxCSRFactDirect( const fax_csr_t      *graphA,
                               const plu_order_t *order,
                               fax_csr_t            *graphL );
Sparse_long faxCSRFactILUk( const fax_csr_t      *graphA,
                             const plu_order_t *order,
                             Sparse_long          level,
                             fax_csr_t            *graphL );

void faxCSRAmalgamate( int             ilu,
                       double          rat_cblk,
                       double          rat_blas,
                       fax_csr_t      *graphL,
                       plu_order_t *order);

/**
 * @}
 */
#endif /* _fax_csr_h_ */

#ifndef _extracblk_h_
#define _extracblk_h_

/**
 * @brief Extra symbol cblk structure
 */
typedef struct extracblk_s {
    Sparse_long   cblknbr; /**< Number of cblk allocated                          */
    Sparse_long   addcblk; /**< Number of cblk created                            */
    Sparse_long   addblok; /**< Number of blok created                            */
    Sparse_long   addblof; /**< Number of blok created due to facing cblk splited */
    Sparse_long  *sptcblk; /**< Index for splitted cblk in the cblktab            */
    Sparse_long  *sptcbnb; /**< Number of splitted cblk for a cblk                */
    Sparse_long   curcblk; /**< Cursor for cblktab                                */
    Sparse_long   sizcblk; /**< Size of allocated cblktab                         */
    symbol_cblk_t *cblktab; /**< Array of column blocks [+1,based]                 */
} ExtraCblk_t;

void extraCblkInit ( Sparse_long        cblknbr,
                     ExtraCblk_t        *extracblk );
void extraCblkExit ( ExtraCblk_t        *extracblk );
void extraCblkAdd  ( ExtraCblk_t        *extracblk,
                     Sparse_long        fcolnum,
                     Sparse_long        lcolnum,
                     int8_t              selevtx );
void extraCblkMerge( const ExtraCblk_t  *extracblk,
                     symbol_matrix_t    *newsymb,
                     Cand              **candtab   );

#endif /* _extracblk_h_ */

#ifndef _graph_h_
#define _graph_h_

/**
 * @brief Graph structure.
 *
 * This structure describes the adjacency graph of a sparse matrix.
 */
struct plu_graph_s {
    Sparse_long  gN;       /**< Global number of vertices in compressed graph    */
    Sparse_long  n;        /**< Number of local vertices in compressed graph     */
    Sparse_long *colptr;   /**< List of indirections to rows for each vertex     */
    Sparse_long *rows;     /**< List of edges for each vertex                    */
    Sparse_long *loc2glob; /**< Corresponding numbering from local to global     */
    Sparse_long *glob2loc; /**< Corresponding numbering from global to local     */
    Sparse_long  dof;      /**< Degre of freedom to move to uncompressed graph   */
    Sparse_long *dofs;     /**< Array of the first column of each element in the
                                 expanded matrix [+1,based]                       */
};

/**
 * @name Graph basic subroutines
 * @{
 */
int  graphPrepare(       plu_data_t   *plu_data,
                   const sparse_csc      *spm,
                         plu_graph_t  **graph );
void graphBase   (       plu_graph_t  *graph, int baseval );
void graphExit   (       plu_graph_t  *graph );

/**
 * @}
 * @name Graph IO subroutines
 * @{
 */
void graphLoad( const plu_data_t  *plu_data,
                plu_graph_t       *graph );
void graphSave( plu_data_t        *plu_data,
                const plu_graph_t *graph );

/**
 * @}
 * @name Graph manipulation subroutines
 * @{
 */
int  graphCopy      ( plu_graph_t       *graphdst,
                      const plu_graph_t *graphsrc );
void graphSort      ( plu_graph_t       *graph );
void graphNoDiag    ( plu_graph_t       *graph );
int  graphSymmetrize(       Sparse_long    n,
                      const Sparse_long   *ia,
                      const Sparse_long   *ja,
                      const Sparse_long   *loc2glob,
                            plu_graph_t *newgraph );

int  graphIsolate   (       Sparse_long    n,
                      const Sparse_long   *colptr,
                      const Sparse_long   *rows,
                            Sparse_long    isolate_n,
                            Sparse_long   *isolate_list,
                            Sparse_long   **new_colptr,
                            Sparse_long   **new_rows,
                            Sparse_long   **new_perm,
                            Sparse_long   **new_invp );

int  graphApplyPerm ( const plu_graph_t *graphA,
                      const Sparse_long   *perm,
                            plu_graph_t *graphPA );

int graphIsolateRange( const plu_graph_t *graphIn,
                       const plu_order_t *order,
		             plu_graph_t *graphOut,
                             Sparse_long    fnode,
                             Sparse_long    lnode,
                             Sparse_long    distance );
void graphComputeProjection( const plu_graph_t *graph,
                             const int            *vertlvl,
                                   plu_order_t *order,
                             const plu_graph_t *subgraph,
                                   plu_order_t *suborder,
                                   Sparse_long    fnode,
                                   Sparse_long    lnode,
                                   Sparse_long    sn_level,
                                   Sparse_long    distance,
                                   Sparse_long    maxdepth,
                                   Sparse_long    maxwidth,
                                   Sparse_long   *depthsze );

Sparse_long graphIsolateConnectedComponents( const plu_graph_t *graph,
                                              Sparse_long         *comp_vtx,
                                              Sparse_long         *comp_sze );

int graphComputeKway( const plu_graph_t *graph,
                      plu_order_t       *order,
                      Sparse_long         *peritab,
                      Sparse_long         *comp_nbr,
                      Sparse_long         *comp_sze,
                      Sparse_long         *comp_vtx,
                      Sparse_long          comp_id,
                      Sparse_long          nbpart );

#endif /* _graph_h_ */

#ifndef _plu_order_h_
#define _plu_order_h_

struct etree_s;
typedef struct etree_s EliminTree;

/**
 * @brief Order structure.
 *
 * This structure stores the permutation (and inverse permutation) associated to the ordering.
 * It also stores the partitioning tree and the set of supernodes.
 */
typedef struct plu_order_s {
    Sparse_long  baseval;   /**< base value used for numbering       */
    Sparse_long  vertnbr;   /**< Number of vertices                  */
    Sparse_long  cblknbr;   /**< Number of column blocks             */
    Sparse_long *permtab;   /**< Permutation array of size vertnbr [based]           */
    Sparse_long *peritab;   /**< Inverse permutation array of size vertnbr [based]   */
    Sparse_long *rangtab;   /**< Supernode array of size cblknbr+1 [based,+1]        */
    Sparse_long *treetab;   /**< Partitioning tree of size cblknbr+1 [based]         */
    int8_t       *selevtx;   /**< Selected vertices for low-rank clustering of size cblknbr */
    Sparse_long  sndenbr;   /**< The number of original supernodes before clustering */
    Sparse_long *sndetab;   /**< Original supernode array of size sndenbr [based,+1] */
} plu_order_t;

/**
 * @name Order basic subroutines
 * @{
 */
int  pluOrderInit  (       plu_order_t * const ordeptr,
                              Sparse_long           baseval,
                              Sparse_long           vertnbr,
                              Sparse_long           cblknbr,
                              Sparse_long   * const perm,
                              Sparse_long   * const invp,
                              Sparse_long   * const rang,
                              Sparse_long   * const tree );
int  pluOrderAlloc (       plu_order_t * const ordeptr,
                              Sparse_long           vertnbr,
                              Sparse_long           cblknbr );
int  pluOrderAllocId(      plu_order_t * const ordeptr,
                              Sparse_long           vertnbr );
void pluOrderExit  (       plu_order_t * const ordeptr );
void pluOrderBase  (       plu_order_t * const ordeptr,
                              Sparse_long           baseval );
int  pluOrderCheck ( const plu_order_t * const ordeptr );
void pluOrderExpand(       plu_order_t * const ordeptr,
                              sparse_csc     * const spm);
int  pluOrderCopy  (       plu_order_t * const ordedst,
                        const plu_order_t * const ordesrc );

const plu_order_t *pluOrderGet( const plu_data_t * const plu_data );

/**
 * @name Order IO subroutines
 */
int  pluOrderLoad( const plu_data_t *plu_data,       plu_order_t *ordeptr );
int  pluOrderSave(       plu_data_t *plu_data, const plu_order_t *ordeptr );

/**
 * @name Order compute subroutines
 */
int  pluOrderComputeMetis(    plu_data_t *plu_data, plu_graph_t *graph );

int  pluOrderGrid( plu_order_t **myorder, Sparse_long nx,
                      Sparse_long ny, Sparse_long nz );

/**
 * @name Order manipulation subroutines
 */
void pluOrderFindSupernodes( const plu_graph_t *graph,
                                plu_order_t * const ordeptr );

int  pluOrderAmalgamate( int             verbose,
                            int             ilu,
                            int             levelk,
                            int             rat_cblk,
                            int             rat_blas,
                            plu_graph_t *graph,
                            plu_order_t *orderptr);

int  pluOrderApplyLevelOrder( plu_order_t *ordeptr,
                                 Sparse_long    level_tasks2d,
                                 Sparse_long    width_tasks2d );

int  pluOrderAddIsolate( plu_order_t     *ordeptr,
                            Sparse_long        new_n,
                            const Sparse_long *perm );

void orderDraw( plu_data_t *plu_data,
                Sparse_long   min_cblk );

Sparse_long
orderSupernodes( const plu_graph_t *graph,
                 plu_order_t       *order,
                 EliminTree           *etree,
                 Sparse_long         *iparm );

#endif /* _plu_order_h_ */

#ifndef _perf_h_
#define _perf_h_

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#define PERF_MODEL "AMD 6180  MKL"

/**GEMM**/
#define GEMM_A  2.429169e-10
#define GEMM_B  2.724804e-10
#define GEMM_C  1.328900e-09
#define GEMM_D  1.148989e-07
#define GEMM_E -2.704179e-10
#define GEMM_F  1.216278e-06
#define PERF_GEMM(i,j,k) (GEMM_A*(double)(i)*(double)(j)*(double)(k)+GEMM_B*(double)(i)*(double)(j)+GEMM_C*(double)(j)*(double)(k)+GEMM_D*(double)(i)+GEMM_E*(double)(j)+GEMM_F)


/**GEAM**/
#define GEAM_A   1.358111e-09
#define GEAM_B  -4.416379e-09
#define GEAM_C   2.270780e-08
#define GEAM_D  -3.335563e-07
#define PERF_GEAM(i,j)   (GEAM_A*(double)(i)*(double)(j)+GEAM_B*(double)(i)+GEAM_C*(double)(j)+GEAM_D)

/**TRSM (Works only for right case) **/
#define TRSM_A 2.626177e-10
#define TRSM_B 3.976198e-08
#define TRSM_C 3.255168e-06
#define PERF_TRSM( i, j )   (TRSM_A*(double)(i)*(double)(i)*(double)(j)+TRSM_B*(double)(i)+TRSM_C)

/**POTRF**/
#define POTRF_A  2.439599e-11
#define POTRF_B  1.707006e-08
#define POTRF_C -1.469893e-07
#define POTRF_D  4.071507e-07
#define PERF_POTRF(i) (POTRF_A*(double)(i)*(double)(i)*(double)(i)+POTRF_B*(double)(i)*(double)(i)+POTRF_C*(double)(i)+POTRF_D)

/**PPF**/
#define PPF_A  2.439599e-11
#define PPF_B  1.707006e-08
#define PPF_C -1.469893e-07
#define PPF_D  4.071507e-07
#define PERF_SYTRF(i) (PPF_A*(double)(i)*(double)(i)*(double)(i)+PPF_B*(double)(i)*(double)(i)+PPF_C*(double)(i)+PPF_D)

/**SCAL**/
#define SCAL_A 4.371793e-10
#define SCAL_B 2.052399e-07
#define PERF_SCAL(i) (SCAL_A*(double)(i)+SCAL_B)

/**COPY**/
#define COPY_A 9.177969e-10
#define COPY_B 2.266129e-07
#define PERF_COPY(i) (COPY_A*(double)(i)+COPY_B)

/**AXPY**/
#define AXPY_A 4.620143e-10
#define AXPY_B 2.101008e-07
#define PERF_AXPY(i) (AXPY_A*(double)(i)+AXPY_B)

/**GEMV**/
#define GEMV_A  6.192657e-10
#define GEMV_B -2.884799e-09
#define GEMV_C  7.594831e-10
#define GEMV_D  3.575035e-07
#define PERF_GEMV(i,j)   (GEMV_A*(double)(i)*(double)(j)+GEMV_B*(double)(i)+GEMV_C*(double)(j)+GEMV_D)

/**TRSV**/
#define TRSV_A 3.224536e-10
#define TRSV_B 1.709178e-08
#define TRSV_C 1.947268e-07
#define PERF_TRSV(i) (TRSV_A*(double)(i)*(double)(i)+TRSV_B*(double)(i)+TRSV_C)

/* en octets ...
   TIME : entre threads */

/* en octets ...
   CLUSTER : entre noeuds */

/* en octets ...
   SHARED : entre MPI shared */

/* old version compatibility
#define TIME_BANDWIDTH    1.5e-9
#define TIME_STARTUP      5.2e-6
#define CLUSTER_BANDWIDTH 5.9e-10
#define CLUSTER_STARTUP   3.9e-6
   end old                  */

#define TIME_BANDWIDTH_1    0.0
#define TIME_STARTUP_1      1e-8
#define SHARED_BANDWIDTH_1  1.0e-10
#define SHARED_STARTUP_1    0.2e-6
#define CLUSTER_BANDWIDTH_1 3.0e-10
#define CLUSTER_STARTUP_1   3.0e-6

#define TIME_BANDWIDTH_2    0.0
#define TIME_STARTUP_2      1e-8
#define SHARED_BANDWIDTH_2  3.0e-10
#define SHARED_STARTUP_2    0.4e-6
#define CLUSTER_BANDWIDTH_2 6.0e-10
#define CLUSTER_STARTUP_2   6.0e-6

#define TIME_BANDWIDTH_4    0.0
#define TIME_STARTUP_4      1e-8
#define SHARED_BANDWIDTH_4  6.0e-10
#define SHARED_STARTUP_4    0.8e-6
#define CLUSTER_BANDWIDTH_4 9.0e-10
#define CLUSTER_STARTUP_4   9.0e-6

#define TIME_BANDWIDTH_8    0.0
#define TIME_STARTUP_8      1e-8
#define SHARED_BANDWIDTH_8  6.0e-10
#define SHARED_STARTUP_8    0.8e-6
#define CLUSTER_BANDWIDTH_8 9.0e-10
#define CLUSTER_STARTUP_8   0.0e-6

#define PENALTY_STARTUP     0.0
#define PENALTY_BANDWIDTH   0.0

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#endif /* _perf_h_ */

