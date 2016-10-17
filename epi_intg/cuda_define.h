/*
 * cuda_define.h
 *
 *  Created on: 2016/10/14
 *      Author: yasu7890v
 */

#ifndef CUDA_DEFINE_H_
#define CUDA_DEFINE_H_
#include "define.h"

#ifdef USE_DOUBLE
typedef double4 real4;
#else
typedef float4 real4;
#endif

#define G_MAX_CONNECT_CELL_NUM (407) //lowest 8n+7 s.t >= 400



#endif /* CUDA_DEFINE_H_ */
