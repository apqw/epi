#pragma once
#include <memory>
#include "define.h"
#define DIST_SQ(x,y,z) ((x)*(x)+(y)*(y)+(z)*(z))

inline double p_diff_x(double x1, double x2)
{
	using namespace cont;
	double diff = x1 - x2;
	if (diff > 0.5*LX)return diff - LX;
	if (diff <= -0.5*LX)return diff + LX;
	return diff;
}

inline double p_diff_y(double y1, double y2)
{
	using namespace cont;
	double diff = y1 - y2;
	if (diff > 0.5*LY)return diff - LY;
	if (diff <= -0.5*LY)return diff + LY;
	return diff;
}

inline double p_dist_sq(double x1, double y1, double z1, double x2, double y2, double z2)
{
	double diffx = p_diff_x(x1, x2);
	double diffy = p_diff_y(y1, y2);
	double diffz = z1 - z2;
	return diffx*diffx + diffy*diffy + diffz*diffz;
}

inline double min0(double a) {
	return a > 0 ? a : 0;
}




double min0(double a);

extern int __lat_x[cont::NX * 3];
extern int __lat_y[cont::NY * 3];
extern int __lat_z[cont::NZ * 3];
extern int per_x_prev_idx[cont::NX];
extern int per_x_next_idx[cont::NX];
extern int per_y_prev_idx[cont::NY];
extern int per_y_next_idx[cont::NY];

void init_precalc_lat();
constexpr int(&precalc_lat_x())[cont::NX * 3]{
	return __lat_x;
}
constexpr int(&precalc_lat_y())[cont::NY * 3]{
	return __lat_y;
}
constexpr int(&precalc_lat_z())[cont::NZ * 3]{
	return __lat_z;
}

void init_precalc_per();

constexpr int(&precalc_per_next_x())[cont::NX]{
	return per_x_next_idx;
}

constexpr int(&precalc_per_prev_x())[cont::NX]{
	return per_x_prev_idx;
}
constexpr int(&precalc_per_next_y())[cont::NY]{
	return per_y_next_idx;
}

constexpr int(&precalc_per_prev_y())[cont::NY]{
	return per_y_prev_idx;
}