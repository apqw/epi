#pragma once
#include <memory>
#include "define.h"

#define DIST_SQ(x,y,z) ((x)*(x)+(y)*(y)+(z)*(z))

double p_diff_x(double x1, double x2);
double p_diff_y(double y1, double y2);
double p_dist_sq(double x1, double y1, double z1, double x2, double y2, double z2);
double p_cell_dist_sq(const Cell* c1, const Cell* c2);
double min0(double a);
bool no_double_count(Cell* c1, Cell* c2);

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