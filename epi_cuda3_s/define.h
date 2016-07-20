#pragma once
#include <cuda.h>
#include <cuda_runtime.h>
#include <thrust/device_ptr.h>
#define __CCONC(x,y) x##y
#define CCONC(x,y) __CCONC(x,y)
#define EVAL(x) x
#define USE_FLOAT

#define USE_INACCURATE_GJ
#define DBG

#ifdef USE_FLOAT

typedef float real;
typedef float4 real4;
typedef float3 real3;
#define CDEF(x) (CCONC(x,F))
#define M_PI_R (3.141592654f)
#define make_real4(...) make_float4(__VA_ARGS__)
#define make_real3(...) make_float3(__VA_ARGS__)
#define sincosr(...) sincosf(__VA_ARGS__)
#define rsqrtr(...) rsqrtf(__VA_ARGS__)
#define sqrtr(...) sqrtf(__VA_ARGS__)
#define R_FMT "%f"
#else
typedef double real;
typedef double4 real4;
typedef double3 real3;
#define CDEF(x) (x)
#define M_PI_R (3.141592653589793238463)
#define make_real4(...) make_double4(__VA_ARGS__)
#define make_real3(...) make_double3(__VA_ARGS__)
#define sincosr(...) sincos(__VA_ARGS__)
#define rsqrtr(...) rsqrt(__VA_ARGS__)
#define sqrtr(...) sqrt(__VA_ARGS__)
#define R_FMT "%lf"
#endif


typedef unsigned int CELL_STATE;
typedef real4 CellPos;
typedef float4 CellPosFP32;
typedef int CellIndex;
typedef float* FloatArr;
typedef real* RealArr;
typedef int* CellMap1;
typedef float* CellMap2;

#define MAX_CELL_NUM (65536u)
#define MAX_CONNECT_CELL_NUM (256u)

#define LX_val 50.0
#define LY_val 50.0
#define LZ_val 100.0

#define LX CDEF(LX_val)
#define LY CDEF(LY_val)
#define LZ CDEF(LZ_val)

#define LXf (CCONC(LX_val,f))
#define LYf (CCONC(LY_val,f))
#define LZf (CCONC(LZ_val,f))

#define NX (100)
#define NY (100)
#define NZ (200)

#define dx (LX/NX)
#define dy (LY/NY)
#define dz (LZ/NZ)

#define inv_dx (NX/LX)
#define inv_dy (NY/LY)
#define inv_dz (NZ/LZ)

#define R_max CDEF(1.4)
#define R_der CDEF(1.4)
#define R_memb CDEF(1.0)

#define LJ_THRESH CDEF(1.2)

#define NMX (100)
#define NMY (100)

#define DT_Cell CDEF(0.01)
#define DT_Ca CDEF(0.01)

#define COMPRESS_FACTOR (4u)

#define THRESH_SP CDEF(3.0)

#define ca2p_init CDEF(0.122)

#define ex_inert_init CDEF(0.97)

#define IP3_init CDEF(0.0)

#define gj_init CDEF(0.99)

#define ATP_init CDEF(0.0)

#define ext_stim_init CDEF(0.0)

template<typename T>
using devPtr = thrust::device_ptr<T>;

template<typename T>
using devRef = thrust::device_reference<T>;

