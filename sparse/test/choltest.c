/* 并行计算平台-Cholesky分解组件 */

#include <sys/time.h>
#include "Sparse.h"
#include "tpsm.h"
#include<math.h>
/**
 * @brief 根据矩阵文件所在根目录以及矩阵名称得到完整路径
 * 
 * @param file_path 
 * @param matrix_name 
 * @param matrix 
 * @return int 
 */
double check_error (sparse_csc *A, sparse_factor *L, sparse_common *cc)
{
    int i, j;
    double one [2] = {1,0}, minusone [2] = {-1,0}, zero [2] = {0,0} ;
    dense_array *X, *B, *X_sol; // 原始X，AX=B， 求解得到X_sol

    int n = A->ncol;
    // n = ones (n,1) 一个稠密的向量: 给稠密矩阵分配空间并置为0 
    X = SparseCore_zeros (n, 1, A->xtype, cc) ;
    B = SparseCore_zeros (n, 1, A->xtype, cc) ;

    // 初始化X = [0, 1, 2,..., n]
    double *x ;
    x = (double *) (X->x) ;
    for (i = 0 ; i < n ; i++)
    {
        x [i] = i ;
    }
    // 计算B = A*X
    SparseCore_sdmult (A, 0, one, zero, X, B, cc) ; 
    
    // 求解X_sol
    X_sol = SparseChol_solve (SPARSE_A, L, B, cc) ;
    
    // 计算两个X的误差
    double diff_norm = 0.0;
    double xx = 0.0;
    double *x_check = (double *) (X_sol->x);
    for (j = 0; j < X->nrow; j++)
    {
        xx =  x_check [j] - j;
        diff_norm += xx * xx;
    }
    diff_norm = sqrt(diff_norm)/X->nrow;

    SparseCore_free_dense (&X, cc) ;
    SparseCore_free_dense (&B, cc) ;
    SparseCore_free_dense (&X_sol, cc) ;
    return diff_norm;
}

int main (int argc, char* argv[])
{
    int cycleNum = 1;
    double timeStart, timeEnd;
    struct timeval tv;
    double one [2] = {1,0}, m1 [2] = {-1,0} ;
    double *Bx ;
    double switch_value;
    sparse_csc *A , *L_sparse = NULL ;
    dense_array *x, *b, *r;
    sparse_factor *L ;
    sparse_common c ;
    etree_info *EtreeInfo = (etree_info*) malloc(sizeof(etree_info));
    SparseCore_start (&c) ;			            /* start HNUCHOL */

    char *fmatrix = argv[1];
    // 提取松弛因子
    // double nzrelax = atof(argv[2]);
    double gamma = argc >= 3 ? atof(argv[2]) : RELAXED_GAMMA;    // 默认1.0
    double relax_cdensity = argc >= 4 ? atof(argv[3]) : RELAXED_CDENSITY;// 默认1.25

    printf("%s\n", fmatrix);
    FILE *fp;
    fp = fopen(fmatrix,"r");
    A = SparseCore_read_sparse (fp, &c) ;	            /* read in a matrix */
    // SparseCore_print_sparse (A, "A", &c) ;		    /* print the matrix */

    /* 结果输出到文件中 */
    char *result_file = (char *)"./ARMcholesky.txt";
    FILE* fresult;
    fresult = fopen(result_file,"a+");
    fprintf(fresult, "%s\t", fmatrix);

    int i;
    if (A == NULL || A->stype == 0)		    /* A must be symmetric */
    {
        printf("Matrix must be symmetric !\n");
	    SparseCore_free_sparse (&A, &c) ;  
	    SparseCore_finish (&c) ;
	return (0) ;
    }

    // 开线程池
    TPSM_init(500, 2000, 1000, TPSM_NODE_AFFINITY );
    c.final_asis = 1;                   // 0 1
    c.supernodal = SPARSE_SUPERNODAL;         //  SPARSE_AUTO SPARSE_SUPERNODAL SPARSE_SIMPLICIAL
    Relaxfactor_setting (A->nrow, SparseCore_nnz (A, &c), RELAX_FOR_CHOL, &c);
    // if(nzrelax!=0.0){
    //     c.nzrelax = nzrelax;
    // }
    // c.nzrelax = nzrelax;
    c.gamma = gamma;
    c.relax_cdensity = relax_cdensity;
    //------------analyze------------
    gettimeofday(&tv, NULL);
    timeStart = tv.tv_sec + tv.tv_usec / 1000000.0;

    L = SparseChol_analyze_p2 (TRUE, A, NULL, NULL, 0, &c, EtreeInfo); /* analyze */

    gettimeofday(&tv, NULL);
    switch_value = c.fl/c.lnz;
    printf("switch_value = %f\n GAMMA=%f, RELAXED=%f\n", switch_value, c.gamma, c.relax_cdensity);
    fprintf(fresult, "%lg\t", c.fl);
    // fprintf(fresult, "%f\t", switch_value);
    // fprintf(fresult, "%d\t", EtreeInfo->nfsuper);
    // fprintf(fresult, "%d\t", EtreeInfo->singlenum);
    // fprintf(fresult, "%f\t", EtreeInfo->avgcols);
    // fprintf(fresult, "%d\t", EtreeInfo->max_depth);
    // fprintf(fresult, "%d\t", EtreeInfo->leafnum);

    timeEnd = tv.tv_sec + tv.tv_usec / 1000000.0;
    double analyze_time = timeEnd - timeStart;
    printf("\n       analyze run time: %f\n", analyze_time);
    // fprintf(fresult, "%lf\t", analyze_time);
    //------------factorize------------
    gettimeofday(&tv, NULL);
    timeStart = tv.tv_sec + tv.tv_usec / 1000000.0;
    
    SparseChol_factorize ( A, L, &c) ;          /* factorize */
    
    gettimeofday(&tv, NULL);
    timeEnd = tv.tv_sec + tv.tv_usec / 1000000.0;
    double factorize_time = timeEnd - timeStart;
    printf("     factorize run time: %f\n", factorize_time);
    printf("SparseCholSuper run time: %f\n\n", analyze_time + factorize_time);
    // printf("SparseCholAuto run time: %f\n\n", analyze_time + factorize_time);
    // fprintf(fresult, "%lf\t", factorize_time);
    // fprintf(fresult, "%.2f\t", c.nzrelax);
    fprintf( fresult, "%.4f\t%.4f\t", c.gamma, c.relax_cdensity);
    fprintf(fresult, "%lf\t%lf\t%lf\t", analyze_time, factorize_time, analyze_time + factorize_time);
    // 等待全部完成,销毁线程池
    TPSM_destroy(TPSM_SHUTDOWN_GENTLY);
    //------------solve------------
    double res = 0.;
    res = check_error(A, L, &c);
    printf ("res = %8.1e \n ", res) ;
    fprintf (fresult, "%8.1e\n", res) ;

    free(EtreeInfo);
    SparseCore_free_factor (&L, &c) ;		    /* free matrices */
    SparseCore_free_sparse (&A, &c) ;
    SparseCore_free_sparse (&L_sparse, &c) ;
    SparseCore_finish (&c) ;			    /* finish HNUCHOL */
    fclose(fp);
    fclose(fresult);
    return (0) ;
}
