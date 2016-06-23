#pragma once
#include <memory>
#include <type_traits>
#include "define.h"
#define DIST_SQ(x,y,z) ((x)*(x)+(y)*(y)+(z)*(z))

/**
*  8近傍の平均を計算
*  @param [in] grd 対象となる3次元配列状のオブジェクト(grd[l][m][n]でアクセスできるオブジェクト)
*  @param [in] ix x座標のインデックス
*  @param [in] iy y座標のインデックス
*  @param [in] iz z座標のインデックス
*  @return (計算できた場合)平均の値
*
*  @note テンプレート関数なのでgrdの型、返り値の型は推定される。
*  ここで平均というのは8つ足して1/8を掛けるだけなので、operator+,operator*が定義されていないと計算できない。
*/
template<typename ArrTy>
inline auto grid_avg8(ArrTy&& grd, int ix, int iy, int iz)->typename std::remove_reference<decltype(grd[ix][iy][iz])>::type { //for c++11
    return 0.125*(grd[ix][iy][iz] + grd[ix + 1][iy][iz] + grd[ix][iy + 1][iz]
        + grd[ix][iy][iz + 1] + grd[ix + 1][iy + 1][iz] + grd[ix + 1][iy][iz + 1]
        + grd[ix][iy + 1][iz + 1] + grd[ix + 1][iy + 1][iz + 1]);
}

inline double p_diff_x(const double x1, const double x2)
{
	using namespace cont;
	const double diff = x1 - x2;
	if (diff > 0.5*LX)return diff - LX;
	if (diff <= -0.5*LX)return diff + LX;
	return diff;
}

inline double p_diff_y(const double y1, const double y2)
{
	using namespace cont;
	const double diff = y1 - y2;
	if (diff > 0.5*LY)return diff - LY;
	if (diff <= -0.5*LY)return diff + LY;
	return diff;
}

inline double p_dist_sq(const double x1, const double y1, const double z1, const double x2, const double y2, const double z2)
{
	const double diffx = p_diff_x(x1, x2);
	const double diffy = p_diff_y(y1, y2);
	const double diffz = z1 - z2;
	return diffx*diffx + diffy*diffy + diffz*diffz;
}

inline double min0(const double a) {
	return a > 0 ? a : 0;
}
/*
template<CELL_STATE_MASK First>
constexpr uint_fast16_t __mask_sum()
{
	return First;
}

template<CELL_STATE_MASK First, CELL_STATE_MASK Second,CELL_STATE_MASK...rest>
constexpr uint_fast16_t __mask_sum()
{
	return First | __mask_sum<Second,rest...>();
}

template<CELL_STATE_MASK... cs>
inline bool state_check(CELL_STATE cellstate) {
	return ((1u << cellstate)&__mask_sum<cs...>()) != 0x0000;
}

*/


double min0(double a);

extern int __lat_x[cont::NX * 3];
extern int __lat_y[cont::NY * 3];
extern int __lat_z[cont::NZ * 3];
extern int per_x_prev_idx[cont::NX];
extern int per_x_next_idx[cont::NX];
extern int per_y_prev_idx[cont::NY];
extern int per_y_next_idx[cont::NY];
extern int per_z_prev_idx[cont::NZ];
extern int per_z_next_idx[cont::NZ];

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

constexpr int(&precalc_per_next_z())[cont::NZ]{
	return per_z_next_idx;
}

constexpr int(&precalc_per_prev_z())[cont::NZ]{
	return per_z_prev_idx;
}