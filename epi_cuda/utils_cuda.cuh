#pragma once
#include "global_cuda.h"
__device__ inline float p_diff_x(const float x1, const float x2)
{
	using namespace cont;
	const float diff = x1 - x2;
	if (diff > 0.5f*LX)return diff - LX;
	if (diff <= -0.5f*LX)return diff + LX;
	return diff;
}

__device__ inline float p_diff_y(const float y1, const float y2)
{
	using namespace cont;
	const float diff = y1 - y2;
	if (diff > 0.5f*LY)return diff - LY;
	if (diff <= -0.5f*LY)return diff + LY;
	return diff;
}
__device__ inline float min0(const float a) {
	return a > 0.0f ? a : 0.0f;
}

__device__ inline float get_radius(unsigned int state){
	using namespace cont;
	return state == MEMB ? R_memb : R_max;
}

__device__ inline float p_cell_dist_sq_d(cell_pos_set a, cell_pos_set b){
	float diffx = p_diff_x(a.x, b.x);
	float diffy = p_diff_y(a.y, b.y);
	float diffz = a.z - b.z;
	return diffx*diffx + diffy*diffy + diffz*diffz;
}