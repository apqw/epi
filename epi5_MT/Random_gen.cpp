#include "Random_gen.h"
#include <algorithm>
#include <array>
Random_gen::Random_gen() :dist(0.0, 1.0) {
	init();
}

void Random_gen::init() {
	std::random_device rng;
	std::array<uint32_t, 16> rng_seed;
	std::generate(rng_seed.begin(), rng_seed.end(), std::ref(rng));
	std::seed_seq rng_seed_seq(rng_seed.begin(), rng_seed.end());
	mt.seed(rng_seed_seq);
}

double Random_gen::gen_rand_real() {
	return dist(mt);
}

double genrand_real()
{
	static thread_local Random_gen rng;
	return rng.gen_rand_real();
}
