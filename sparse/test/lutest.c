/* 并行计算平台-LU分解组件 */

#include "Sparse.h"

int main (int argc, char* argv[])
{
    char *fmatrix;
    int LUmethod;
    Sparse_long iparm[IPARM_SIZE];            /*< Integer in/out parameters for plu                */
    double dparm[DPARM_SIZE];               /*< Floating in/out parameters for plu               */

    // Initialize parameters to default values
    pluInitParam(iparm, dparm);
    // 读取命令行参数
    pluGetOptions(argc, argv, &LUmethod, iparm, dparm, NULL, &fmatrix);
    printf("%s\n", fmatrix);

    /* 结果输出到文件中 */
    char *result_file = (char *)"./ARM_LU.txt";
    FILE* fresult;
    fresult = fopen(result_file,"a+");
    fprintf(fresult, "%s\t", fmatrix);

    FILE *fp;
    fp = fopen(fmatrix,"r");
    // 如果读取矩阵不存在，返回错误
    if(fp == NULL) {
        fprintf(fresult, "No exist\n");
        printf("Error: %s file is not exist!\n", fmatrix);
        free(fmatrix);
        fclose(fresult);
        return 0;
    }
    /* 根据矩阵类型、行列数、非零元数量等信息自适应选择LU分解方法 */
    adaptiveLU(LUmethod, fp, fmatrix, fresult, iparm, dparm);
    free(fmatrix);
    fclose(fp);
    fclose(fresult);
    return (0) ;
}
