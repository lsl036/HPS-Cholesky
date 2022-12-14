2020年11月22日更新记录:

    1. 新增了Makefile.option 所有相关的 makefile 变量可以在其中修改。

    2. 替换了2.2版本的线程池，使得并行计算的时间更加稳定。

    3. 对于所有的稀疏矩阵分解操作，我们调整了超节点合并的松弛因子，并加入了函数能参考矩阵的稀疏度来调整合并的条件。经过部分数据的测试，目前可以提升一些稀疏矩阵的分解性能。

    4. 对于稀疏矩阵 QR 分解，调整了稠密矩阵分解时的分块参数：FCHUNK, SMALL等，对于QR分解提升了一些性能。

    5. 对于Cholesky分解，我们新增了 METIS 填充简化排序方法，并在/lib 文件夹下加入了 libmetis.a 静态库，这是由 metis 函数包编译而来，该函数包是Apache协议的，符合使用条件。
        在 Cholesky 分解的过程中，会首先使用 AMD 重排序方法分析，然后根据AMD得出的 flops 和 nnz 的比值，当大于 1500 的时候会开启 METIS 方法再做符号分析，能够进一步的缩小填充元素数目，在增加符号分解时间的同时也能减少数值分解的时间。 这个阈值 1500 是选用了一部分数据后的测试结果，在这个分界点条件下，计算的总时间要比单纯使用 AMD 计算的总时间要更快。 新增的方法并不是必须要开启的（开关这一方法的编译参数还没有添加），使用 METIS 方法可以提高一部分数据集的Cholesky分解的效率。

2020年11月30日更新记录:

    1. 添加了CMake编译方式，编译好的库在build目录下。

    2. 为区分两个版本的稀疏LU库，现将之前HNUSparse库下的稀疏LU重命名为olu，意为采用openmp方式并行的LU。将现在新研究的稀疏LU命名为plu，意为采用pthread方式并行的LU。目前plu已合并至HNUSparse库，但尚未研究出自动调用olu或plu的合适阈值。

2020年12月20日更新记录:
    1. 整合了plu头文件。
    2. 找出服务器上不存在的矩阵
    3. 整合OLU和PLU数据结构，统一读命令行函数（参数：读矩阵、线程数、LU分解方法）
    4. 设置读矩阵头部函数，来确定采用何种LU分解算法（统一读矩阵函数）
    6. 设置LU阈值
    7. 编辑LU文档
    8. 补充LU注释

2020年12月24日更新记录:
    1. 修改了Cmake中出现的问题。调整了稠密分解中的连接问题。
    2. 整合了QR、LQ进行测试的最新版本，可以达到与测试报告上的性能。