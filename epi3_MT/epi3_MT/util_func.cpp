#include "util_func.h"
#include "Cell.h"
#define _USE_MATH_DEFINES
#include <math.h>

Random_gen::Random_gen():dist(0.0,1.0){
    init();
}

void Random_gen::init(){
    std::random_device rng;
    std::array<uint32_t, 16> rng_seed;
    std::generate(rng_seed.begin(), rng_seed.end(), std::ref(rng));
    std::seed_seq rng_seed_seq(rng_seed.begin(), rng_seed.end());
    mt.seed(rng_seed_seq);
}

double Random_gen::gen_rand_real(){
return dist(mt);
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
    static thread_local Random_gen rng;
    return rng.gen_rand_real();
}

double p_diff_sc_x(double v1, double v2) {
	using namespace cont;
	double diff = v1 - v2;
	if (diff > 0.5*LX)return diff - LX;
	if (diff <= -0.5*LX)return diff + LX;
	return diff;

}

double p_diff_sc_y(double v1, double v2) {
	using namespace cont;
	double diff = v1 - v2;
	if (diff > 0.5*LY)return diff - LY;
	if (diff <= -0.5*LY)return diff + LY;
	return diff;

}

Vec3<DV<double>> p_diff_v3(Vec3<DV<double>> V1, Vec3<DV<double>> V2) {
	return Vec3<DV<double>>({
		p_diff_sc_x((double)V1[0],(double)V2[0]),
		p_diff_sc_y((double)V1[1],(double)V2[1]),
		V1[2] - V2[2]
	});
}

double normSqV3(Vec3<DV<double>> V1) {
	return V1[0]() * V1[0]() + V1[1]() * V1[1]() + V1[2]() * V1[2]();
}

double distSqV3(Vec3<DV<double>> V1, Vec3<DV<double>> V2) {
	double diffx = p_diff_sc_x(V1[0](), V2[0]());
	double diffy = p_diff_sc_y(V1[1](), V2[1]());
	double diffz = V1[2]() - V2[2]();

	return diffx*diffx + diffy*diffy + diffz*diffz;
}

double cellDistSq(Cell* c1, Cell *c2) {
	return distSqV3(c1->pos, c2->pos);
}

bool is_near(Cell* c1, Cell* c2) {
	return distSqV3(c1->pos, c2->pos) < (c1->radius + c2->radius)*(c1->radius + c2->radius);
}

bool is_near_delta(Cell* c1, Cell* c2, double delta) {
	return sqrt(distSqV3(c1->pos, c2->pos)) < (c1->radius + c2->radius) - delta;
}

double ljmain(Cell* c1, Cell* c2) {
	return _ljmain(c1->radius(), c2->radius(), cellDistSq(c1, c2));
}

double ljmain_der_near(Cell* c1, Cell* c2) {
	return _ljmain_der_near(c1->radius(), c2->radius(), sqrt(cellDistSq(c1, c2)));
}

double ljmain_der_far(Cell* c1, Cell* c2) {
	return _ljmain_der_far(c1->radius(), c2->radius(), sqrt(cellDistSq(c1, c2)));
}

double adhesion(Cell* c1, Cell* c2, double sprconst) {
	return _adhe(sqrt(cellDistSq(c1, c2)), c1->radius() + c2->radius(), sprconst);
}

bool paired_with_fix(Cell* c1) {
	return c1->pair != nullptr&&c1->pair->state() == FIX;
}

Vec3<double> calc_dermis_normal(Cell * me, Cell * dermis)
{
	Vec3<double> nv = me->pos - dermis->pos;
	double norm = sqrt(vec_sum(nv*nv));
	return nv / norm;
}

Vec3<double> div_direction(Cell* me, Cell* dermis) {
	auto normal = calc_dermis_normal(me, dermis);
	double rand_theta, rand_phi, cr1, sr1, cr2, sr2;
	Vec3<double> randv = { 0,0,0 };
	Vec3<double> out = { 0,0,0 };
    double sumSq;
	do {
		rand_theta = M_PI*genrand_real();
		cr1 = cos(rand_theta);
		sr1 = sin(rand_theta);
		rand_phi = 2 * M_PI*genrand_real();
		cr2 = cos(rand_phi);
		sr2 = sin(rand_phi);
		randv[0] = sr1*cr2;
		randv[1] = sr1*sr2;
		randv[2] = cr1;

		out = cross(normal, randv);

	} while ((sumSq = vec_sum(out*out)) < 1.0e-14);
	return out / sqrt(sumSq);
}

uint32_t BitSeparateFor3D(uint32_t n)
{
	uint32_t s = n;
	s = (s | s << 8) & 0x0000f00f;
	s = (s | s << 4) & 0x000c30c3;
	s = (s | s << 2) & 0x00249249;
	return s;
}


uint32_t Get3DMortonOrder(uint32_t x, uint32_t y, uint32_t z)
{
	return BitSeparateFor3D(x) | BitSeparateFor3D(y) << 1 | BitSeparateFor3D(z) << 2;
}
