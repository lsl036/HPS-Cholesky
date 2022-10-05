#ifndef _solver_h_
#define _solver_h_

struct blendctrl_s;
typedef struct blendctrl_s BlendCtrl;

struct simuctrl_s;
typedef struct simuctrl_s SimuCtrl;

#include "plu_dlrcores.h"

/**
 * @name Cblk properties
 *  The type and structure definitions.
 *  Define the mask for the cblks in the cblktype field:
 *   - 1st bit: The cblk is a fake local cblk corresponding to a fanin that will be sent to someone else
 *   - 2nd bit: The cblk is stored in a 2D layout fashion as in a tiled matrix, otherwise the standard 1D lapack layout is used
 *   - 3rd bit: The cblk generates 2D granularity tasks, instead of a single 1D tasks that perform factorization, solves and updates
 *   - 4th bit: The cblk is compressed in Low-Rank (implies CBLK_LAYOUT_2D), otherwise it is stored in dense
 *   - 5th bit: The cblk is part of the Schur complement if set
 */
#define CBLK_FANIN      (1 << 0)
#define CBLK_LAYOUT_2D  (1 << 1)
#define CBLK_TASKS_2D   (1 << 2)
#define CBLK_COMPRESSED (1 << 3)
#define CBLK_IN_SCHUR   (1 << 4)
#define CBLK_IN_LAST    (1 << 5)

/*
 * The type and structure definitions.
 */
#define COMP_1D                     0
#define DIAG                        1
#define E1                          2
#define E2                          3
#define DRUNK                       4

/**
 * @brief The task structure for the numerical factorization
 */
typedef struct task_s {
    Sparse_long          taskid;  /**< COMP_1D DIAG E1 E2                                        */
    Sparse_long          prionum; /**< Priority value for the factorization                      */
    Sparse_long          cblknum; /**< Attached column block                                     */
    Sparse_long          bloknum; /**< Attached block                                            */
    Sparse_long volatile ftgtcnt; /**< Number of fan-in targets                                  */
    Sparse_long volatile ctrbcnt; /**< Total number of contributions                             */
    Sparse_long          indnum;  /**< For E2 (COMP_1D), index of ftgt (>0) else if local = -taskdest
                                        For DIAG and E1, index of btag (>0) if there is a
                                        local one it must be the first of the chain of local E1   */
#if defined(PLU_DYNSCHED)
    Sparse_long          threadid;/**< Index of the bubble which contains the task               */
#endif
} Task;

/**
 * @brief Fan-in target information fields
 * @warning The number of fields must be odd for memory alignment purpose.
 */
typedef enum {
    FTGT_CTRBNBR = 0,           /**< Number of contributions            */
    FTGT_CTRBCNT,               /**< Number of contributions remaining  */
    FTGT_PROCDST,               /**< Destination for fanintarget        */
    FTGT_TASKDST,               /**< Task  destination                  */
    FTGT_BLOKDST,               /**< Block destination (->COMP_1D)      */
    FTGT_PRIONUM,               /**< Fanintarget priority               */
    FTGT_FCOLNUM,               /**< Fanintarget first column           */
    FTGT_LCOLNUM,               /**< Fanintarget last column            */
    FTGT_FROWNUM,               /**< Fanintarget first row              */
    FTGT_LROWNUM,               /**< Fanintarget last row               */
    FTGT_MAXINFO
} solver_ftgt_e;

/**
 * @brief Fan-in target structure for data exchange
 */
typedef struct solver_ftgt_s {
    Sparse_long   infotab[FTGT_MAXINFO]; /**< Fan-in target header holding all information enumerated in solver_ftgt_e */
    void          *coeftab;               /**< Fan-in target coeficient array                                           */
} solver_ftgt_t;

/**
 * @brief Solver block structure.
 */
typedef struct solver_blok_s {
    void        *handler[2]; /**< Runtime data handler                     */
    Sparse_long lcblknm;    /**< Local column block                       */
    Sparse_long fcblknm;    /**< Facing column block                      */
    Sparse_long frownum;    /**< First row index                          */
    Sparse_long lrownum;    /**< Last row index (inclusive)               */
    Sparse_long coefind;    /**< Index in coeftab                         */
    Sparse_long browind;    /**< Index in browtab                         */
    int8_t       inlast;     /**< Index of the block among last separator (2), coupling with last separator (1) or other blocks (0) */

    /* LR structures */
    plu_lrblock_t *LRblock; /**< Store the blok (L/U) in LR format. Allocated for the cblk. */
} SolverBlok;

/**
 * @brief Solver column block structure.
 */
typedef struct solver_cblk_s  {
    plu_atomic_lock_t lock;       /**< Lock to protect computation on the cblk */
    volatile uint32_t    ctrbcnt;    /**< Number of contribution to receive       */
    int8_t               cblktype;   /**< Type of cblk                            */
    Sparse_long         fcolnum;    /**< First column index                      */
    Sparse_long         lcolnum;    /**< Last column index (inclusive)           */
    SolverBlok          *fblokptr;   /**< First block in column (diagonal)        */
    Sparse_long         stride;     /**< Column block stride                     */
    Sparse_long         lcolidx;    /**< Local first column index to the location in the rhs vector       */
    Sparse_long         brownum;    /**< First block in row facing the diagonal block in browtab, 0-based */
    Sparse_long         brown2d;    /**< First 2D-block in row facing the diagonal block in browtab, 0-based */
    Sparse_long         gcblknum;   /**< Global column block index               */
    Sparse_long         sndeidx;    /**< Index of the original supernode the cblk belongs to */
    void                *lcoeftab;   /**< Coefficients access vector              */
    void                *ucoeftab;   /**< Coefficients access vector              */
    void                *handler[2]; /**< Runtime data handler                    */
    Sparse_long         selevtx;    /**< Index to identify selected cblk for which intra-separator contributions are not compressed */
    Sparse_long         threadid;   /**< Rank of the accessing thread            */
} SolverCblk;

/**
 * @brief Solver column block structure.
 *
 * This structure stores all the numerical information about the factorization,
 * as well as the structure of the problem. Only local information to each
 * process is stored in this structure.
 *
 */
struct solver_matrix_s {
    int restore; /**< Flag to indicate if it is require to restore data with
                      solverBackupRestore: 0: No need, 1:After solve,
                      2:After Factorization */
    Sparse_long            baseval;       /**< Base value for numberings                         */
    Sparse_long            nodenbr;       /**< Number of nodes before dof extension              */
    Sparse_long            coefnbr;       /**< Number of coefficients (node after dof extension) */
    Sparse_long            gcblknbr;      /**< Global number of column blocks                    */
    Sparse_long            cblknbr;       /**< Number of column blocks                   */
    Sparse_long            cblkmax1d;     /**< Rank of the last cblk not beeing enabled for 2D computations */
    Sparse_long            cblkmin2d;     /**< Rank of the first cblk beeing enabled for 2D computations        */
    Sparse_long            cblkmaxblk;    /**< Maximum number of blocks per cblk         */
    Sparse_long            cblkschur;     /**< Index of the first local cblk in Schur    */
    Sparse_long            nb2dcblk;      /**< Number of 2D cblks                        */
    Sparse_long            nb2dblok;      /**< Number of 2D blocks                       */
    Sparse_long            bloknbr;       /**< Number of blocks                          */
    Sparse_long            brownbr;       /**< Size of the browtab array                 */
    SolverCblk   * restrict cblktab;       /**< Array of solver column blocks             */
    SolverBlok   * restrict bloktab;       /**< Array of solver blocks                    */
    Sparse_long * restrict browtab;       /**< Array of blocks                           */

    plu_lr_t             lowrank;       /**< Low-rank parameters                       */
    plu_factotype_t      factotype;     /**< General or symmetric factorization?       */
    double                  diagthreshold; /**< Diagonal threshold for pivoting           */
    volatile int32_t        nbpivots;      /**< Number of pivots during the factorization */

    Sparse_long              ftgtnbr;              /*+ Number of fanintargets                    +*/
    Sparse_long              ftgtcnt;              /*+ Number of fanintargets to receive         +*/
    solver_ftgt_t * restrict  ftgttab;              /*+ Fanintarget access vector                 +*/

    Sparse_long              offdmax;              /*+ Maximum size of the off-diagonal blocks for hetrf/sytrf temporary buffers +*/
    Sparse_long              gemmmax;              /*+ Maximum size of the GEMM update for 1d GEMM computations                  +*/
    Sparse_long              blokmax;              /*+ Maximum size of 2D blocks                 +*/
    Sparse_long              nbftmax;              /*+ Maximum block number in ftgt              +*/
    Sparse_long              arftmax;              /*+ Maximum block area in ftgt                +*/

    Sparse_long              clustnum;             /*+ current processor number                  +*/
    Sparse_long              clustnbr;             /*+ number of processors                      +*/
    Sparse_long              procnbr;              /*+ Number of physical processor used         +*/
    Sparse_long              thrdnbr;              /*+ Number of local computation threads       +*/
    Sparse_long              bublnbr;              /*+ Number of local computation threads       +*/
    /* BubbleTree   * restrict   btree;                /\*+ Bubbles tree                              +*\/ */

    Sparse_long              indnbr;
    Sparse_long * restrict   indtab;
    Task         * restrict   tasktab;              /*+ Task access vector                        +*/
    Sparse_long              tasknbr;              /*+ Number of Tasks                           +*/
    Sparse_long **           ttsktab;              /*+ Task access vector by thread              +*/
    Sparse_long *            ttsknbr;              /*+ Number of tasks by thread                 +*/
    plu_queue_t **         computeQueue;         /*+ Queue of task to compute by thread        +*/

    Sparse_long *            proc2clust;           /*+ proc -> cluster                           +*/
    Sparse_long              gridldim;             /*+ Dimensions of the virtual processors      +*/
    Sparse_long              gridcdim;             /*+ grid if dense end block                   +*/

    Sparse_long             *selevtx;              /*+ Array to identify which cblk are pre-selected +*/

};

/**
 * @brief     Compute the number of columns in a column block.
 * @param[in] cblk
 *            The pointer to the column block.
 * @return    The number of columns in the cblk.
 */
static inline Sparse_long
cblk_colnbr( const SolverCblk *cblk )
{
    return cblk->lcolnum - cblk->fcolnum + 1;
}

/**
 * @brief     Compute the number of blocks in a column block.
 * @param[in] cblk
 *            The pointer to the column block.
 * @return    The number of blocks in the cblk including the diagonal block.
 */
static inline Sparse_long
cblk_bloknbr( const SolverCblk *cblk )
{
    return (cblk+1)->fblokptr - cblk->fblokptr + 1;
}

/**
 * @brief     Compute the number of rows of a block.
 * @param[in] blok
 *            The pointer to the block.
 * @return    The number of rows in the block.
 */
static inline Sparse_long
blok_rownbr( const SolverBlok *blok )
{
    return blok->lrownum - blok->frownum + 1;
}

/**
 * @brief     Compute the number of rows of a column block.
 * @param[in] cblk
 *            The pointer to the column block.
 * @return    The number of rows in the column block.
 */
static inline Sparse_long
cblk_rownbr( const SolverCblk *cblk )
{
    Sparse_long rownbr = 0;
    SolverBlok * blok;
    for (blok = cblk->fblokptr; blok < cblk[1].fblokptr; blok++) {
        rownbr += blok_rownbr(blok);
    }
    return rownbr;
}

/**
 * @brief    Task stealing method.
 *
 * @param[inout] solvmtx
 *            The pointer to the solverMatrix.
 * @param[in] rank
 *            Rank of the computeQueue.
 * @param[inout] dest
 *            Rank of the stolen queue.
 * @param[in] nbthreads
 *            Total amount of threads in the node.
 */
static inline Sparse_long
stealQueue( SolverMatrix *solvmtx,
            int           rank,
            int          *dest,
            int           nbthreads )
{
    int rk = *dest;
    plu_queue_t *stoleQueue;
    Sparse_long    cblknum = -1;
    while( rk != rank )
    {
        assert( solvmtx->computeQueue[ rk ] );
        stoleQueue = solvmtx->computeQueue[ rk ];
        if( (cblknum = pqueuePop(stoleQueue)) != -1 ){
            *dest = rk;
            return cblknum;
        }
        rk = (rk + 1)%nbthreads;
    }
    return cblknum;
}

/**
 * @brief Check if a block is included inside another one.
 *
 * Indicate if a blok is included inside an other block.
 * i.e. indicate if the row range of the first block is included in the
 * one of the second.
 *
 * @param[in] blok  The block that is tested for inclusion.
 * @param[in] fblok The block that is suppose to include the first one.
 *
 * @retval true   if the first block is     included in the second one.
 * @retval false  if the first block is not included in the second one.
 */
static inline int
is_block_inside_fblock( const SolverBlok *blok,
                        const SolverBlok *fblok )
{
#  if defined(NAPA_SOPALIN)
    return (((blok->frownum >= fblok->frownum) &&
             (blok->lrownum <= fblok->lrownum)) ||
            ((blok->frownum <= fblok->frownum) &&
             (blok->lrownum >= fblok->lrownum)) ||
            ((blok->frownum <= fblok->frownum) &&
             (blok->lrownum >= fblok->frownum)) ||
            ((blok->frownum <= fblok->lrownum) &&
             (blok->lrownum >= fblok->lrownum)));
#  else
    return ((blok->frownum >= fblok->frownum) &&
            (blok->lrownum <= fblok->lrownum));
#  endif /* defined(NAPA_SOPALIN) */
}

void solverInit( SolverMatrix *solvmtx );
void solverExit( SolverMatrix *solvmtx );

int  solverMatrixGen( Sparse_long           clustnum,
                      SolverMatrix          *solvmtx,
                      const symbol_matrix_t *symbmtx,
                      const plu_order_t  *ordeptr,
                      const SimuCtrl        *simuctl,
                      const BlendCtrl       *ctrl );

int  solverLoad( SolverMatrix       *solvptr,
                 FILE               *stream );
int  solverSave( const SolverMatrix *solvptr,
                 FILE               *stream );

void          solverRealloc( SolverMatrix       *solvptr);
SolverMatrix *solverCopy   ( const SolverMatrix *solvptr,
                             int                 flttype );

int           solverCheck     ( const SolverMatrix *solvmtx );
int           solverDraw      ( const SolverMatrix *solvptr,
                                FILE               *stream,
                                int                 verbose,
                                const char         *directory );
void          solverPrintStats( const SolverMatrix *solvptr );

/*
 * Solver backup
 */
struct SolverBackup_s;
typedef struct SolverBackup_s SolverBackup_t;

SolverBackup_t *solverBackupInit   ( const SolverMatrix *solvmtx                          );
int             solverBackupRestore(       SolverMatrix *solvmtx, const SolverBackup_t *b );
void            solverBackupExit   (                                    SolverBackup_t *b );

#endif /* _solver_h_ */
