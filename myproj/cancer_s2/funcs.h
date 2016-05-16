#ifndef FUNCS_H
#define FUNCS_H
#include "define.h"
#include <random>
static bool rng_initalized = false;
static std::mt19937_64 mt;
static std::uniform_real_distribution<double> rng_real_dist(0.0, 1.0);

void genrand_init();
double genrand_real();

double p_diff_sc_x(double v1, double v2);

double p_diff_sc_y(double v1, double v2);

double p_dist_sq(double x1,double y1,double z1,double x2,double y2,double z2);

double calc_LJ_force_coef(double distSq,double radius1,double radius2);
#endif // FUNCS_H
