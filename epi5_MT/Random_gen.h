#pragma once
#include <random>
class Random_gen {
	std::mt19937_64 mt;
	std::uniform_real_distribution<double> dist;
public:
	void init();
	Random_gen();
	double gen_rand_real();
};

double genrand_real();