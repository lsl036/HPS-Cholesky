'''
Description: 
Author: Shengle Lin
Date: 2022-04-23 12:53:36
LastEditors: Shengle Lin
LastEditTime: 2022-09-14 09:50:11
'''
import os
import sys
import numpy as np

# -----------------------------------
# --     执行程序choltest_nnz      --
# -----------------------------------
def runtest( func_name, file_name, option):
    cmd1 = 'export OMP_NUM_THREADS=32'
    cmd2 = "./%s %s %s" %(func_name, file_name, option)
    cmd = cmd1 + "&&" + cmd2
    val = os.system(cmd)


def runalldata (file_name):
    print('--- Runtest Begin ---')
    # file_name = 'sample.txt'
    func_name = 'build/test/sparse_choltest_nnz'
    total_data = open(file_name)
    for line in total_data.readlines():
        formLine = line.strip().split(' ')
        dataname = formLine[0]
        border = float( formLine[1])
        GAMMA = np.linspace(-1*border,border,num=5)
        # GAMMA = np.linspace(-1*border,border,num=3) 
        RELAXED = np.linspace(1.05,1.5,num=10)
        for gamma in GAMMA:
            for relaxed in RELAXED:
                option = '%f %f' %(gamma, relaxed)
                runtest(func_name, dataname, option)
    print('--- Runtest Finished ---')

# 记录最优性能的 alpha, beta 选择
def dataprocess (filename):
    print('--- Process Data Begin ---')
    Fp = open(filename)
    last_line = Fp.readline()
    formlastLine = last_line.strip().split('\t')
    min_data = [ formlastLine[0], formlastLine[-3], formlastLine[-2], formlastLine[-1]]
    final_result=[]
    for line in Fp.readlines():
        currentLine = line.strip().split('\t')
        # 如果是同一个数据集的结果, 才进行比较
        if currentLine[0] == min_data[0]:
            if float(currentLine[-1]) < float(min_data[-1]): # 这一行的时间更短
                min_data[-3] = currentLine[-3]
                min_data[-2] = currentLine[-2]
                min_data[-1] = currentLine[-1]
        else:   # 新的数据集开始了，写之后存新的
            # 保存当前数据的最小值结果
            final_result.append(min_data)
            # 写入新的数据集信息
            min_data = []
            min_data = [ currentLine[0], currentLine[-3], currentLine[-2], currentLine[-1]]
    final_result.append(min_data)
    dst_data = open('proceed.txt','a+')

    for output in final_result:
        dst_data.write("{} {} {} {}".format(output[0],output[1],output[2],output[3])+'\n')
    print('--- Process Data Finished ---')

# 用于筛选出未完成计算的数据集
def supplement_data(file_all, file_finished):
    result = []
    Fp_all   = open(file_all)
    Fp_ready = open(file_finished)
    lines_all = Fp_all.readlines()
    lines_ready = Fp_ready.readlines()
    # result = list (set(lines_all).difference(set(Fp_ready)))
    for line_t in lines_all:
        if line_t not in lines_ready:
            result.append(line_t)
    dst = open("supplementdata0516.txt",'a+')
    dst.writelines(result)
    # for output in result:
    #     dst.write()

# 原始生成的数据集存在问题，保存了全部的非零元， clean_data 用于生成仅保留三角部分的数据集
def clean_data(file_list):
    print('--- Clean Data Begin ---')
    File_listp = open(file_list)
    for line in File_listp.readlines():
        line = line.strip()
        file_name = "/home/public/GEN_SYM/" + line +".mtx"
        dst_name = "/home/public/GEN_CLEAN/" + line +".mtx"
        # 打开原始的数据集
        Fp = open(file_name,"r+")
        Dstp = open(dst_name,'a+')
        #  读写第一行的 banner
        bannerline = Fp.readline()
        Dstp.write("{}".format(bannerline))
        # 读写第二行 head
        headLine = Fp.readline().strip().split(' ')
        m = int(headLine[0])
        n = int(headLine[1])
        nnz = int(headLine[2])
        file_nnz = int((nnz + m)/2)
        print(file_nnz)
        Dstp.write("{} {} {}".format(m,n,file_nnz)+'\n')
        for dataline in Fp.readlines():
            data = dataline.strip().split(' ')
            if int(data[0])>= int(data[1]):
                Dstp.write("{}".format(dataline))
        Dstp.close()
    print('--- Clean Data Finish ---')

if __name__ == '__main__':

    # clean_data('filelist.txt')
    # 脚本， 大规模测试数据集
    # runalldata('resdata0516.txt')

    # 找每个数据集对应的最优解
    # dataprocess('ARMcholesky.txt')

    # 找出还未测试的数据集
    # supplement_data('allgendata.txt', 'finished_dataset.txt')

    # print('--- Runtest begin ---')
    # funcname = 'build/test/sparse_choltest'
    

    # 批量进行数据测试
    # totalfile = open('dataset.txt')
    # for line in totalfile.readlines():
    #     line = line.strip('\n')
    #     file_name = '/home/public/DATA_TEST/'+ line +'/'+ line + '.mtx'
    #     # print(file_name)
    #     GAMMA = 1.0
    #     RELAXED = np.linspace(80,100,num=21) 
    #     for relaxed in RELAXED:
    #         option = '%f %f' %(GAMMA, relaxed)
    #         runtest(funcname, file_name, option)

    # GAMMA = np.linspace(1,1.5,num=11)
    # RELAXED = np.linspace(0.8,1.5,num=71)

    # 单独数据测试
    # filename = '/home/public/DATA_TEST/thermal2/thermal2.mtx'
    # GAMMA = np.linspace(1,5,num=51)
    # RELAXED = np.linspace(1.0,1.5,num=21)  
    # for gamma in GAMMA:
    # for relaxed in RELAXED:
    #     option = '%f %f' %(0.363653884, relaxed)
    #     runtest(funcname, filename, option)
    # val = os.system('ls -al')
    # print ("val = %d" %(val))
