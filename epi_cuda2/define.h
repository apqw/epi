/*
 * define.h
 *
 *  Created on: 2016/07/07
 *      Author: yasu7890v
 */

#ifndef DEFINE_H_
#define DEFINE_H_
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#ifdef _WIN32
#define RESTRICT 
#else
#define RESTRICT 
#endif

typedef unsigned int CELL_STATE;
typedef float4 CellPos;
typedef int CellIndex;
typedef float* FloatArr;

#define MAX_CELL_NUM (65536u)
#define MAX_CONNECT_CELL_NUM (200u)

#define LX (100.0f)
#define LY (50.0f)
#define LZ (100.0f)

#define R_max (1.4f)
#define R_der (1.4f)
#define R_memb (1.0f)

#define LJ_THRESH (1.2f)

#define NMX (300)
#define NMY (150)

#define DT_Cell (0.01f)

#define COMPRESS_FACTOR (6u)

#define THRESH_SP (3.0f)

#define ca2p_init (0.122f)

#define ex_inert_init (0.97f)

#define IP3_init (0.0f)

#define gj_init (0.99f)

#define ATP_init (0.0f)

#define ext_stim_init (0.0f)

enum __CST:unsigned int {
    ALIVE = 0u,
    DEAD = 1u,
    DISA = 2u,
    UNUSED = 3u, //
    FIX = 4u,
    BLANK = 5u,
    DER = 6u,
    MUSUME = 7u,
    AIR = 8u,
    MEMB = 9u
};

#define SYSTEM 0
#define WHOLE 0
#define BASAL 1
#define STOCHASTIC 1

#define MALIGNANT (0)
#define M_PI_F (3.141592654f)

#define THRESH_DEAD (22.0f)

#endif /* DEFINE_H_ */
