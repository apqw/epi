#include "funcs.h"
#include <cassert>
double ret_x(double x1, double x2)
/* choose one of LX-periodic values of x2 such that
|x1-x2| takes minimum
*/
{
	double min, tmp;
	min = fabs(x2 - x1);
	tmp = x2;
	if (min > fabs(x2 - LX - x1)) {
		min = fabs(x2 - LX - x1);
		tmp = x2 - LX;
	}
	if (min > fabs(x2 + LX - x1)) {
		min = fabs(x2 + LX - x1);
		tmp = x2 + LX;
	}
	return tmp;
}
double ret_y(double x1, double x2)
{
	double min, tmp;
	min = fabs(x2 - x1);
	tmp = x2;
	if (min > fabs(x2 - LY - x1)) {
		min = fabs(x2 - LY - x1);
		tmp = x2 - LY;
	}
	if (min > fabs(x2 + LY - x1)) {
		min = fabs(x2 + LY - x1);
		tmp = x2 + LY;
	}
	return tmp;
}
void genrand_init() {
	std::random_device rng;
	std::array<uint32_t, 16> rng_seed;
	std::generate(rng_seed.begin(), rng_seed.end(), std::ref(rng));
	std::seed_seq rng_seed_seq(rng_seed.begin(), rng_seed.end());
	mt.seed(rng_seed_seq);
	rng_initalized = true;
}
double genrand_real() {
	assert(rng_initalized);
	return rng_real_dist(mt);
}

double bend_force_sqr(
	const Arr1D<double>& n,
	const Arr1D<double>& m,
	const Arr1D<double>& ipn,
	const Arr1D<double>& ipm,
	const Arr1D<double>& dn,
	const Arr1D<double>& dm,
	int j, int jr, int jl, int jll, int ju, int jb, int jbb)
{
	return KBEND * (
		-(1.0 - ipn[j])*(n[jr] - ipn[j] * n[j]) / dn[j]
		+ (1.0 - ipn[jl])*(-(n[jl] - ipn[jl] * n[j]) / dn[j]
			+ (n[j] - ipn[jl] * n[jl]) / dn[jl])
		+ (1.0 - ipn[jll])*(n[jll] - ipn[jll] * n[jl]) / dn[jl]

		- (1.0 - ipm[j])*(m[ju] - ipm[j] * m[j]) / dm[j]
		+ (1.0 - ipm[jb])*(-(m[jb] - ipm[jb] * m[j]) / dm[j]
			+ (m[j] - ipm[jb] * m[jb]) / dm[jb])
		+ (1.0 - ipm[jbb])*(m[jbb] - ipm[jbb] * m[jb]) / dm[jb]);
}