/*
 * utils_cuda.h
 *
 *  Created on: 2016/07/07
 *      Author: yasu7890v
 */
#ifndef UTILS_CUDA_H_
#define UTILS_CUDA_H_
#include "define.h"
#include <string>
#include <sstream>

enum DIRECTION{
	DIR_U,DIR_L,DIR_B,DIR_R,NONE
};
template<CELL_STATE cs>
__host__ __device__ inline real get_radius(){
	return R_max;
}
template<>
__host__ __device__ inline real get_radius<MEMB>(){
	return R_memb;
}
template<>
__host__ __device__ inline real get_radius<DER>(){
	return R_der;
}

__host__ __device__ inline real get_radius(CELL_STATE cs){
	return cs == MEMB ? R_memb : cs == DER ? R_der : R_max;
}



__host__ __device__ inline real p_diff_x(const real x1,const real x2){
	const real diff= x1-x2;
	if(diff>0.5f*LX)return diff-LX;
	if(diff<=-0.5f*LX)return diff+LX;
	return diff;
}

__host__ __device__ inline float p_diff_x_f(const float x1, const float x2){
	const float diff = x1 - x2;
	if (diff>0.5f*LXf)return diff - LXf;
	if (diff <= -0.5f*LXf)return diff + LXf;
	return diff;
}

__host__ __device__ inline real p_diff_y(const real y1,const real y2){
	const real diff = y1-y2;
	if(diff>0.5f*LY)return diff-LY;
	if(diff<=-0.5f*LY)return diff+LY;
	return diff;
}

__host__ __device__ inline float p_diff_y_f(const float y1, const float y2){
	const float diff = y1 - y2;
	if (diff>0.5f*LYf)return diff - LYf;
	if (diff <= -0.5f*LYf)return diff + LYf;
	return diff;
}

__host__ __device__ inline real min0(const real a){
	return a>0.0?a:0.0;
}

__host__ __device__ inline float min0f(const float a){
	return a>0.0f ? a : 0.0f;
}

__host__ __device__ inline real p_cell_dist_sq(const CellPos a,const CellPos b){
	const real diffx = p_diff_x(a.x,b.x);
	const real diffy = p_diff_y(a.y,b.y);
	const real diffz = a.z-b.z;
	return diffx*diffx + diffy*diffy + diffz*diffz;
}

__host__ __device__ inline float p_cell_dist_sq_f(const CellPosFP32 a, const CellPosFP32 b){
	const float diffx = p_diff_x(a.x, b.x);
	const float diffy = p_diff_y(a.y, b.y);
	const float diffz = a.z - b.z;
	return diffx*diffx + diffy*diffy + diffz*diffz;
}
template<unsigned X,unsigned Y,unsigned Z>
__host__ __device__ inline int midx(int i,int j,int k){
	return i*Y*Z+j*Z+k;
}


template<DIRECTION d>
__device__ __host__ inline CellIndex get_adj_memb_idx(CellIndex idx){
	return 0;
}//dummy


template<>
__device__ __host__ inline CellIndex get_adj_memb_idx<DIR_U>(CellIndex my_memb_idx){
	const int kk= my_memb_idx / NMX;
	return ((kk + 1) % NMY)*NMX + my_memb_idx%NMX;
}

template<>
__device__ __host__ inline CellIndex get_adj_memb_idx<DIR_L>(CellIndex my_memb_idx){
	const int jj = my_memb_idx%NMX;
	return ((int)(my_memb_idx / NMX))*NMX + (jj - 1+NMX) % NMX;
}

template<>
__device__ __host__ inline CellIndex get_adj_memb_idx<DIR_B>(CellIndex my_memb_idx){
	const int kk = my_memb_idx / NMX;
	return ((kk - 1+NMY) % NMY)*NMX + my_memb_idx%NMX;
}

template<>
__device__ __host__ inline CellIndex get_adj_memb_idx<DIR_R>(CellIndex my_memb_idx){
	const int jj = my_memb_idx%NMX;
	return ((int)(my_memb_idx / NMX))*NMX + (jj + 1) % NMX;
}


__device__ __host__ inline real4 operator+(const real4 a, const real4 b){
	return make_real4(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
}

__device__ __host__ inline real4 operator-(const real4 a, const real4 b){
	return make_real4(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w);
}

__device__ __host__ inline real4 operator*(const real4 a, const real4 b){
	return make_real4(a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w);
}

__device__ __host__ inline real4 operator*(const real c, const real4 b){
	return make_real4(c * b.x, c * b.y, c * b.z, c * b.w);
}

/*
__device__ __host__ inline real sum4(const real4 a){
	return a.x + a.y + a.z + a.w;
}
*/
__device__ __host__ inline real sum3(const real4 a){
	return a.x + a.y + a.z;
}
/*
__device__ __host__ inline real dot4(const real4 a, const real4 b){
	return sum4(a*b);
}
*/
__device__ __host__ inline real dot3(const real4 a, const real4 b){
	return sum3(a*b);
}

__device__ __host__ inline real4 p_diff4(const real4 a, const real4 b){
	return make_real4(p_diff_x(a.x, b.x), p_diff_y(a.y, b.y), a.z - b.z, 0.0f);
}

__device__ __host__ inline real4 p_fms4(const real c, const real4 a, const real4 b){
	return c*p_diff4(a, b);
}
__device__ __host__ inline real interp(real start,real end,real alpha){
	return start*(1.0f-alpha)+end*alpha;
}
std::string int_to_string(int number);

template<size_t X,size_t Y,size_t Z>
__device__ __host__ inline real grid_avg8(const real* grd, int ix, int iy, int iz){
	return 0.125f*(grd[midx<X, Y, Z>(ix, iy, iz)] + grd[midx<X, Y, Z>(ix + 1, iy, iz)] + grd[midx<X, Y, Z>(ix, iy+1, iz)]
		+ grd[midx<X, Y, Z>(ix, iy, iz + 1)] + grd[midx<X, Y, Z>(ix + 1, iy + 1, iz)] + grd[midx<X, Y, Z>(ix+1, iy, iz+1)]
		+ grd[midx<X, Y, Z>(ix, iy + 1, iz + 1)] + grd[midx<X, Y, Z>(ix+1, iy+1, iz+1)]);
}

#if __CUDACC_VER_MAJOR__ < 8 && !defined(USE_FLOAT)
__device__ double atomicAdd(double* address, double val);
#endif


#endif /* UTILS_CUDA_H_ */
