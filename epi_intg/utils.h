#pragma once
#include <string>
#include "util/vec/Vec.h"
#include "define.h"
#define DIST_SQ(x,y,z) ((x)*(x)+(y)*(y)+(z)*(z))

std::string operator"" _s(const char* ,std::size_t);

real p_diff_x(const real x1, const real x2);
real p_diff_y(const real y1, const real y2);
real p_dist_sq(const real x1, const real y1, const real z1, const real x2, const real y2, const real z2);
real min0(const real a);
std::string state_to_str(CELL_STATE);


template<size_t N, typename Intg>
Vec<N, Intg> vpsub(const Vec<N, Intg>& v1, const Vec<N, Intg>& v2) {
    Vec<N, Intg> tmp;
    tmp[0] = p_diff_x(v1[0], v2[0]);
    tmp[1] = p_diff_y(v1[1], v2[1]);
    tmp[2] = v1[2] - v2[2];
    return tmp;
}

template<size_t N,typename Intg>
Vec<N, Intg> vpsub(const Vec<N, Intg>& v1, Vec<N, Intg>&& v2) {
    v2[0] = p_diff_x(v1[0], v2[0]);
    v2[1] = p_diff_y(v1[1], v2[1]);
    v2[2] = v1[2]- v2[2];
    return std::move(v2);
}

template<size_t N, typename Intg>
Vec<N, Intg> vpsub(Vec<N, Intg>&& v1,const Vec<N, Intg>& v2) {
    v1[0] = p_diff_x(v1[0], v2[0]);
    v1[1] = p_diff_y(v1[1], v2[1]);
    v1[2] = v1[2] - v2[2];
    return std::move(v1);
}

template<size_t N, typename Intg>
Vec<N, Intg> vpsub(Vec<N, Intg>&& v1,Vec<N, Intg>&& v2) {
    v1[0] = p_diff_x(v1[0], v2[0]);
    v1[1] = p_diff_y(v1[1], v2[1]);
    v1[2] = v1[2] - v2[2];
    return std::move(v1);
}
