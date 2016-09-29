#pragma once
#include <string>
#include "define.h"
#define DIST_SQ(x,y,z) ((x)*(x)+(y)*(y)+(z)*(z))

std::string operator"" _s(const char* ,std::size_t);

real p_diff_x(const real x1, const real x2);
real p_diff_y(const real y1, const real y2);
real p_dist_sq(const real x1, const real y1, const real z1, const real x2, const real y2, const real z2);
real min0(const real a);