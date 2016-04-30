#pragma once
#include "const.h"
#include <random>
#include <array>
#define SQ(val) ((val)*(val))
static bool rng_initalized = false;
static std::mt19937_64 mt;
static std::uniform_real_distribution<double> rng_real_dist(0.0,1.0);

double ret_x(double x1, double x2);
double ret_y(double x1, double x2);
void genrand_init();
double genrand_real();
double bend_force_sqr(
	const Arr1D<double>& n,
	const Arr1D<double>& m,
	const Arr1D<double>& ipn,
	const Arr1D<double>& ipm,
	const Arr1D<double>& dn,
	const Arr1D<double>& dm,
	int j, int jr, int jl, int jll, int ju, int jb, int jbb);
inline double perio_diff_x(const double x1, const double x2)
{
	double diff = x1 - x2;
	if (diff*diff < 0.25*LX*LX) return diff;
	else {
		return diff + (x1 > x2 ? -LX : LX);
	}
}

inline double perio_diff_y(const double x1, const double x2)
{
	double diff = x1 - x2;
	if (diff*diff < 0.25*LY*LY) return diff;
	else {
		return diff + (x1 > x2 ? -LY : LY);
	}
}

inline double normSq(const double x1, const double x2, const double x3) {
	return x1*x1 + x2*x2 + x3*x3;
}
inline double dist3Sq(const double x1, const double y1, const double z1,
	const double x2, const double y2, const double z2) {
	double diffx = perio_diff_x(x1, x2);
	double diffy = perio_diff_y(y1, y2);
	return normSq(perio_diff_x(x1, x2), perio_diff_y(y1, y2), z1 - z2);
}

inline double adhesion(double distlj, double rad_sum, double spring_const)
{
	double LJ2_m1;

	LJ2_m1 = distlj / rad_sum -1;
	return (LJ2_m1+1 > LJ_THRESH ?
		0.0 :
		-(spring_const / distlj)
		* LJ2_m1*(1 - LJ2_m1*LJ2_m1 / ((LJ_THRESH - 1.0)*(LJ_THRESH - 1.0))));
}