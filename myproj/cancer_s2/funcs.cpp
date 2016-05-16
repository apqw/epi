#include "funcs.h"
void genrand_init() {
    std::random_device rng;
    std::array<uint32_t, 16> rng_seed;
    std::generate(rng_seed.begin(), rng_seed.end(), std::ref(rng));
    std::seed_seq rng_seed_seq(rng_seed.begin(), rng_seed.end());
    mt.seed(rng_seed_seq);
    rng_initalized = true;
}
double genrand_real() {
    //assert(rng_initalized);
    //static thread_local std::mt19937_64 mt;
    //std::uniform_real_distribution<double> distr(0.0,1.0);
    return rng_real_dist(mt);
}

double p_diff_sc_x(double v1, double v2) {
    double diff = v1 - v2;
    if (diff > 0.5*LX)return diff - LX;
    if (diff <= -0.5*LX)return diff + LX;
    return diff;

}

double p_diff_sc_y(double v1, double v2) {
    double diff = v1 - v2;
    if (diff > 0.5*LY)return diff - LY;
    if (diff <= -0.5*LY)return diff + LY;
    return diff;

}

double p_dist_sq(double x1, double y1, double z1, double x2, double y2, double z2){
    double diffx=p_diff_sc_x(x1,x2);
    double diffy=p_diff_sc_y(y1,y2);
    double diffz=z1-z2;
    return diffx*diffx+diffy*diffy+diffz*diffz;
}


/*
 * F=c*dx+c*dy+c*dz ‚Ìc‚ðŒvŽZ
 */
double calc_LJ_force_coef(double distSq,double radius1,double radius2){
    double rad_sum_sq=radius1+radius2;
        rad_sum_sq*=rad_sum_sq;
        double LJ6=rad_sum_sq/distSq;
        LJ6=LJ6*LJ6*LJ6;
return 4.0*LJ_eps*LJ6*(12.0*LJ6-6.0)/distSq;
}
