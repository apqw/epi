#pragma once
#include "define.h"
#include "component.h"
#include "primitive_func.h"
#include <cmath>
#include <random>
class Cell;

static bool rng_initalized = false;
static std::mt19937_64 mt;
static std::uniform_real_distribution<double> rng_real_dist(0.0, 1.0);

void genrand_init();
double genrand_real();

double p_diff_sc_x(double v1, double v2);

double p_diff_sc_y(double v1, double v2);

Vec3<DV<double>> p_diff_v3(Vec3<DV<double>> V1, Vec3<DV<double>> V2);

double normSqV3(Vec3<DV<double>> V1);

double distSqV3(Vec3<DV<double>> V1, Vec3<DV<double>> V2);

double cellDistSq(Cell* c1, Cell *c2);

bool is_near(Cell* c1, Cell* c2);

bool is_near_delta(Cell* c1, Cell* c2, double delta);

double ljmain(Cell* c1, Cell* c2);

double ljmain_der_near(Cell* c1, Cell* c2);

double ljmain_der_far(Cell* c1, Cell* c2);

double adhesion(Cell* c1, Cell* c2, double sprconst);

bool paired_with_fix(Cell* c1);

template<typename T,unsigned N>
T vec_sum(const Vec<T, N>& v) {
	T tmp = v[0];
	for (int i = 1; i < N; i++) {
		tmp += v[i];
	}
	return tmp;
}
template<typename T>
Vec3<T> cross(const Vec3<T>& v1,const Vec3<T>& v2) {
	return Vec3<T>(
	{ 
		v1[1] * v2[2] - v1[2] * v2[1],
		v1[2] * v2[0] - v1[0] * v2[2],
		v1[0] * v2[1] - v1[1] * v2[0] 
	});

}

Vec3<double> calc_dermis_normal(Cell* me, Cell* dermis);
Vec3<double> div_direction(Cell* me, Cell* dermis);

inline CELL_STATE_MASK get_state_mask(CELL_STATE state) {
	return (CELL_STATE_MASK)(1u << (unsigned int)state);
}

inline double grid_avg8(const Arr3D<DV<double>>& grd,int ix,int iy,int iz) {
	return 0.125*(grd[ix][iy][iz]() + grd[ix + 1][iy][iz]() + grd[ix][iy + 1][iz]()
		+ grd[ix][iy][iz + 1]() + grd[ix + 1][iy + 1][iz]() + grd[ix + 1][iy][iz + 1]()
		+ grd[ix][iy + 1][iz + 1]() + grd[ix + 1][iy + 1][iz + 1]());
}

inline double grid_avg8_n(const double grd[cont::NX + 1][cont::NY + 1][cont::NZ+1], int ix, int iy, int iz) {
	return 0.125*(grd[ix][iy][iz] + grd[ix + 1][iy][iz] + grd[ix][iy + 1][iz]
		+ grd[ix][iy][iz + 1] + grd[ix + 1][iy + 1][iz] + grd[ix + 1][iy][iz + 1]
		+ grd[ix][iy + 1][iz + 1] + grd[ix + 1][iy + 1][iz + 1]);
}


template<class T> constexpr T smax(T a) noexcept {
	return a;
}

template<class T> constexpr T smax(T a, T b) noexcept {
	return a > b ? a : b;
}

template<class T, class ... Args>
constexpr T smax(T a, T b, Args ... args) noexcept {
	return smax(a, smax(b, args...));
}

template<class T> constexpr T smin(T a) noexcept {
	return a;
}

template<class T> constexpr T smin(T a, T b) noexcept {
	return a > b ? b : a;
}

template<class T, class ... Args>
constexpr T smin(T a, T b, Args ... args) noexcept {
	return smin(a, smin(b, args...));
}

template<class T>
constexpr T spow(T base, T exp) noexcept {
	//static_assert(exp >= 0, "Exponent must not be negative");
	return exp <= 0 ? 1
		: exp == 1 ? base
		: base * spow(base, exp - 1);
}

uint32_t BitSeparateFor3D(uint32_t n);

uint32_t Get3DMortonOrder(uint32_t, uint32_t, uint32_t);