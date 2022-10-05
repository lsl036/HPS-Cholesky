
#ifndef CHOLMOD_PARTITION_H
#define CHOLMOD_PARTITION_H

#include "SparseCore.h"

#define Int Sparse_long

/* -------------------------------------------------------------------------- */
/* cholmod_metis */
/* -------------------------------------------------------------------------- */

/* Order A, AA', or A(:,f)*A(:,f)' using METIS_NodeND. */

int SparseCore_metis
(
    /* ---- input ---- */
    sparse_csc *A,	/* matrix to order */
    Int *fset,		/* subset of 0:(A->ncol)-1 */
    size_t fsize,	/* size of fset */
    int postorder,	/* if TRUE, follow with etree or coletree postorder */
    /* ---- output --- */
    Int *Perm,		/* size A->nrow, output permutation */
    /* --------------- */
    sparse_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_metis_bisector */
/* -------------------------------------------------------------------------- */

/* Find a set of nodes that bisects the graph of A or AA' (direct interface
 * to METIS_ComputeVertexSeperator). */

Sparse_long SparseCore_metis_bisector	/* returns separator size */
(
    /* ---- input ---- */
    sparse_csc *A,	/* matrix to bisect */
    Int *Anw,		/* size A->nrow, node weights, can be NULL, */
                        /* which means the graph is unweighted. */ 
    Int *Aew,		/* size nz, edge weights (silently ignored). */
                        /* This option was available with METIS 4, but not */
                        /* in METIS 5.  This argument is now unused, but */
                        /* it remains for backward compatibilty, so as not */
                        /* to change the API for cholmod_metis_bisector. */
    /* ---- output --- */
    Int *Partition,	/* size A->nrow */
    /* --------------- */
    sparse_common *Common
) ;

#endif
