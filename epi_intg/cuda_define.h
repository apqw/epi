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




#endif /* CUDA_DEFINE_H_ */
