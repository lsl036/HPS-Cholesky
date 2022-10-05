/*
 * @Description: 
 * @Author: Shengle Lin
 * @Date: 2022-06-06 10:59:06
 * @LastEditors: Shengle Lin
 * @LastEditTime: 2022-09-14 14:39:42
 */
/* 并行计算平台-Cholesky 特征记录组件 */

#include <sys/time.h>
#include "Sparse.h"
#include "tpsm.h"
#include<math.h>
int main (int argc, char* argv[])
{
    int mtype, stype_new, xtype_new;
    sparse_csc *A , *A_read ;
    sparse_factor *L;
    sparse_common c, c_read;
    etree_info *EtreeInfo = (etree_info*) malloc(sizeof(etree_info));
    

    char *fmatrix = argv[1];
    int graph_id = atoi(argv[2]);
    double gamma = argc >= 4 ? atof(argv[3]) : RELAXED_GAMMA;    // 默认1.0
    double relax_cdensity = argc >= 5 ? atof(argv[4]) : RELAXED_CDENSITY;// 默认1.25
    int label = atoi(argv[5]); // label = argv[5]

    printf("---  %s writing  ---\n", fmatrix);

    FILE *fp;
    fp = fopen(fmatrix,"r");
    // A = SparseCore_read_sparse (fp, &c) ;	            /* read in a matrix */
    if(fp == NULL) {
        printf("%s file is not exist!\n", fmatrix);
        return 0;
    }
    #ifdef WRITE_GRAPH
    SparseCore_start (&c) ;			            /* start HNUCHOL */
    A = SparseCore_read_sparse (fp, &c) ;
    fclose(fp);
    long real_nnz = A->nzmax * 2  - A->diag_nz;

    c.final_asis = 1;                   // 0 1
    c.supernodal = SPARSE_SUPERNODAL;
    Relaxfactor_setting (A->nrow, real_nnz, RELAX_FOR_CHOL, &c);

    c.gamma = gamma;
    c.relax_cdensity = relax_cdensity;
    //------------analyze------------

    L = SparseChol_analyze_p2 (TRUE, A, NULL, NULL, 0, &c, EtreeInfo); /* analyze */


    FILE * fresult_Node, *fresult_extinfo;
    char *Node_file = (char *)"./raw/adj_matformat.content";
    char *extinfo_file = (char *)"./raw/matformat_extinfo.txt";
    fresult_Node = fopen(Node_file ,"a+");
    fresult_extinfo = fopen(extinfo_file ,"a+");
    for (size_t k = 0; k < A->nrow; k++)
    {
        /*  Node depth + |R| of node + column width + size of L */
        fprintf(fresult_Node, "%d %ld %ld %ld %ld %ld\n", graph_id, k, (EtreeInfo->node_depth)[k], (EtreeInfo->max_row)[k], (EtreeInfo->nscol)[k], (EtreeInfo->Lsize)[k]);
    }
    fclose(fresult_Node);

    fprintf(fresult_extinfo, "%d %ld %ld %lg %lf %lf %d %d\n", graph_id, A->nrow, real_nnz, (double) real_nnz / (A->nrow*A->nrow), EtreeInfo->switch_value, EtreeInfo->avgcols, EtreeInfo->max_depth, EtreeInfo->leafnum);
    fclose(fresult_extinfo);

    FILE * fresult_Y;
    char *Y_file = (char *)"./raw/matformat_y.txt";
    fresult_Y = fopen(Y_file ,"a+");
    fprintf(fresult_Y, "%d %d\n", graph_id, label);
    fclose(fresult_Y);

    /* Write Edge*/
    fp = fopen(fmatrix,"r");
    SparseCore_start (&c_read) ;			            /* start HNUCHOL */
    A_read = (sparse_csc *) SparseCore_read_matrix (fp, 1, &mtype, &stype_new, &xtype_new, graph_id, EtreeInfo, &c_read);

    SparseCore_free_factor (&L, &c) ;		    /* free matrices */
    SparseCore_free_sparse (&A, &c) ;
    SparseCore_free_sparse (&A_read, &c_read) ;
    SparseCore_finish (&c) ;			    /* finish HNUCHOL */
    SparseCore_finish (&c_read) ;			    /* finish HNUCHOL */
    fclose(fp);
    #else
        printf ("Undefined Write Graph, Please recompile this Library with 'WRITE_GRAPH' definition. \n");
    #endif

    return 0;
}