/**************************************************************
* INCLUDE
**************************************************************/
#include <stdio.h>
#include <pthread.h>
#include <stdlib.h>
#include <unistd.h>
#include "tpsm_base.h"
#include "tpsm_auxiliary.h"
#include "tpsm.h"
/**************************************************************
* DEFINE
**************************************************************/

#define a(i,j)  matrixA[j*lda+i];

pthread_mutex_t lock;


/**************************************************************
* STRUCTURE
**************************************************************/
typedef struct
{
    int * a_row;
	int * b_column;
	int length;
	int ans;
	int father_rank;
	int syntag;
}CALCULATEC_t;


typedef struct
{
    int  a;
    int  b;
    int  ans;
	int father_rank;
}MUL_t;
/**************************************************************
* PRIVATE FUNCTION
**************************************************************/
//get column major matrix element
int get_column_major_matrix_element
(	
	int * const Matrix,
	int const leading_dimention,
	int const row_rank,
	int const col_rank
)
{
	return Matrix[col_rank*leading_dimention+row_rank];
}

//打印矩阵
void printMatrix(int * const Matrix,int const M,int const N, int const ld)
{
	for(int i=0;i<M;++i)
	{
		for(int j=0; j<N; ++j)
		{
			printf("%d ",get_column_major_matrix_element(Matrix,M,i,j));
		}
		printf("\n");
	}
	printf("\n\n");
}


void *mul(void * arg)
{	
	MUL_t *p = (MUL_t*)arg;
	int my_rank = TPSM_get_myrank();	
	int my_nrk = TPSM_Numa_GetNodeRankSelf();
	printf("父线程(%d)\t当前线程(%d)\t函数:%s\t当前节点(%d)\t开始\n",p->father_rank, my_rank, __FUNCTION__,my_nrk);

	p->ans=(p->a)*(p->b);

	printf("父线程(%d)\t当前线程(%d)\t函数:%s\t当前节点(%d)\t完成\n",p->father_rank, my_rank, __FUNCTION__,my_nrk);
}

void *calculateElementC(void * arg)
{	
	int checkcode;
	CALCULATEC_t * calc = (CALCULATEC_t*)arg;
	int const my_rank = TPSM_get_myrank();
	int const my_nrk = TPSM_Numa_GetNodeRankSelf();
	int const syntag = calc->syntag;
	printf("父线程(%d)\t当前线程(%d)\t函数:%s\t当前节点(%d)\t开始\n",calc->father_rank,my_rank,__FUNCTION__,my_nrk);

	MUL_t* mulnode = malloc(calc->length*sizeof(*mulnode));
	for(int i=0;i<calc->length;++i)
	{
		mulnode[i].a = calc->a_row[i];
		mulnode[i].b = calc->b_column[i];
		mulnode[i].ans = 0 ;
		mulnode[i].father_rank = my_rank;
		checkcode = TPSM_addTask(mul,&mulnode[i],syntag);
		TPSM_assert(checkcode, 0);
	}

	printf("父线程(%d)\t当前线程(%d)\t函数:%s\t当前节点(%d)\tTPSM_barrier_tag(2)前\n",calc->father_rank,my_rank,__FUNCTION__,my_nrk);
	TPSM_barrier_tag(syntag);
	printf("父线程(%d)\t当前线程(%d)\t函数:%s\t当前节点(%d)\tTPSM_barrier_tag(2)后\n",calc->father_rank,my_rank,__FUNCTION__,my_nrk);
	for(int i=0;i<calc->length;++i)
	{
		calc->ans += mulnode[i].ans;
	}

	printf("父线程(%d)\t当前线程(%d)\t函数:%s\t当前节点(%d)\t完成\n",calc->father_rank,my_rank,__FUNCTION__,my_nrk);
}




/**************************************************************
* PUBLIC FUNCTION
**************************************************************/

//c = c[0]+c[1]+...+c[NEED_THREADS-1]
int tgemm
(	
	int *  matrixC , 
	int const M, 
	int const N ,
	int const K, 
	int  * 	const matrixA,
	int const lda,
	int  *  const matrixB,
	int const ldb
)
{   
	int checkcode;
	//开辟matrixC元素个数的节点数组
	CALCULATEC_t * calc = malloc(M*N*sizeof(*calc));
	TPSM_assert(calc,1);
	int my_rank = TPSM_get_myrank();
	int my_nrk = TPSM_Numa_GetNodeRankSelf();
	printf("线程(%d)\t 当前节点(%d)\t 执行函数%s\t开始-------\n",my_rank, my_nrk, __FUNCTION__);

	int idx;
	for(int j=0; j<N; ++j)
	{
		for(int i=0; i<M; ++i)
		{	
			idx = j*M+i; 
			
			//存matrixA的i行的个数
			calc[idx].a_row = malloc(K*sizeof(*calc[idx].a_row));
			TPSM_assert(calc[idx].a_row,1);

			//存matrixB的j列的个数
			calc[idx].b_column = malloc(K*sizeof(*calc[idx].b_column));
			TPSM_assert(calc[idx].b_column,1);

			calc[idx].ans = 0;
			calc[idx].length = K;
			calc[idx].father_rank = my_rank;
			calc[idx].syntag = idx+1;
			for(int k=0; k<calc[idx].length;++k)
			{
				(calc[idx].a_row)[k] =  get_column_major_matrix_element(matrixA,lda,i,k);
				(calc[idx].b_column)[k] =  get_column_major_matrix_element(matrixB,ldb,k,j);
				// printf("在%d行%d列的calc[%d]的:\
				// 		a_row[%d]=%d;\
				// 		b_column[%d]=%d\n",\
				// 		i,j,idx,\
				// 		k,(calc[idx].a_row)[k],\
				// 		k,(calc[idx].b_column)[k]);
			}
			//printf("%d\n",idx);
			checkcode = TPSM_addTask(calculateElementC,&calc[idx],0);
			TPSM_assert(checkcode, 0);
			//printf("%d%d\n",idx,idx);
		}
	}

	TPSM_barrier_tag(0);

	for(int i=0; i<M*N; ++i)
	{
		matrixC[i] = calc[i].ans;
	}

	printf("线程(%d)\t 当前节点(%d)\t 执行函数%s\t完成-------\n",my_rank,my_nrk,__FUNCTION__);
  	return 0;
}
/**************************************************************
* MAIN
**************************************************************/
int main(int argc, char **argv)
{   

    int matrixA[16];
    int matrixB[16];
    int matrixC[16];
	int M = 4;
	int N = 4;
	int K = 4;
	int lda = M;
	int ldb = K;
	for(int i=0;i<16;++i)
	{
		matrixA[i] = i+1;
		matrixB[i] = i+21;
		matrixC[i] = 0;
	}

	printf("-------------------------------------\n");
	printf("输入:\n");
	printf("A:\n");
	printMatrix(matrixA,M,K,M);
	printf("B:\n");
	printMatrix(matrixB,K,N,K);
	printf("C:\n");
	printMatrix(matrixC,M,N,M);
	printf("-------------------------------------\n");

    TPSM_init( 128, 100, 100,  TPSM_NODE_AFFINITY);
    
	tgemm(matrixC, M, N, K, matrixA, lda, matrixB, ldb);
    
	printMatrix(matrixC,M,N,M);
    
	TPSM_destroy(TPSM_SHUTDOWN_GENTLY);
	
	int * arr_node0 = TPSM_Numa_SearchNodeSequence(0);
	int * arr_node1 = TPSM_Numa_SearchNodeSequence(1);
	int * arr_node2 = TPSM_Numa_SearchNodeSequence(2);
	int * arr_node3 = TPSM_Numa_SearchNodeSequence(3);
	return 0;
}