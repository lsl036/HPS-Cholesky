/*
 * @Description: 
 * @Author: Shengle Lin
 * @Date: 2022-10-26 13:45:30
 * @LastEditors: Shengle Lin
 * @LastEditTime: 2022-10-26 13:45:56
 */
/********************************************************************************************************************************
 * VERSION: 2.6-Release
 * AUTHOR:  蔡沁耘@湖南大学
 * EMAIL:   hnutsai@hnu.edu.cn
 * DATE:    2020年12月16日
 * FILE:    tpsm_sysinfo.h
 * BRIEF:   TPSM并行框架Release版本(version 2.6),为TPSM_v2.5版本的精炼版本	
 ********************************************************************************************************************************/

#ifndef TPSM_SYSINFO_H
#define TPSM_SYSINFO_H

/******************************************************************************
 * 系统参数 
 *****************************************************************************/

//-------CPU------->
// #define TPSM_CPUSOCKET 		(2) 		//2个cpu_socket (2块cpu)
// #define TPSM_NUMANODES 		(4) 		//4个numa结点
// #define TPSM_SYSCORES 		(128)  		//128核
#include"Numainfo.h"

//-------Cache------->
#define TPSM_L1DOUBLE 		(8192) 		//L1cache可容纳双精度浮点数的个数为：8192(个)
#define TPSM_L1FLOAT 		(16384) 	//L1cache可容纳单精度浮点数的个数为：16384(个)
#define TPSM_CACHELINE 		(128) 		//缓存行大小为：128(Byte字节)

#endif