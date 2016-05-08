#include "primitive_func.h"
double _ljmain(double r1, double r2, double distSq) {
	double rad_sum = r1 + r2;
	double LJ6 = rad_sum*rad_sum / distSq;
	LJ6 = LJ6*LJ6*LJ6;
	return 4.0*cont::eps_m*LJ6*(LJ6 - 1) / distSq;
}

double _ljmain_der_near(double r1, double r2, double dist) {
	double dist2 = dist + cont::delta_R;
	double rad_sum = r1 + r2;
	double LJ6 = rad_sum*rad_sum / (dist2*dist2);
	LJ6 = LJ6*LJ6*LJ6;
	return 4.0*cont::eps_m*LJ6*(LJ6 - 1) / (dist*dist2) + cont::para_ljp2;
}

double _ljmain_der_far(double r1, double r2, double dist) {
	double dist2 = dist + cont::delta_R;
	double rad_sum = r1 + r2;
	double LJ6 = rad_sum*rad_sum / (dist2*dist2);
	LJ6 = LJ6*LJ6*LJ6;
	return 4.0*cont::eps_m*LJ6*(LJ6 - 1) / (dist*dist) + cont::para_ljp2;
}

double _adhe(double distlj, double rad_sum, double spring_const)
{
	using namespace cont;
	double LJ2_m1;

	LJ2_m1 =( distlj / rad_sum) - 1;
	return (LJ2_m1 + 1 > LJ_THRESH ?
		0.0 :
		-(spring_const / distlj)
		* LJ2_m1*(1 - LJ2_m1*LJ2_m1 / ((LJ_THRESH - 1.0)*(LJ_THRESH - 1.0))));
}

double min0(double u) {
	return u > 0 ? u : 0;
}

double fB(double age, double B, bool cornif) {
	using namespace cont;
	return (cornif&&age > THRESH_DEAD - DUR_ALIVE&&age <= THRESH_DEAD + DUR_DEAD ? 1 : 0) - kb*B;
}

double fu(double u, double v, double p, double B)
{
	//  static double result_debug=0.0;
	using namespace cont;
	return kf * (mu0 + mu1 * p / (p + kmu)) * (para_b + para_bb * u / (para_k1 + u)) * v
		- gamma * u / (kg + u) + beta + kbc * B*B * Cout / (Hb + B*B);

}

double fw(double diff, double w)
{
	using namespace cont;
	return (1. - w) + (-1. + tanh((wd - diff) / 0.1)) / 2.; //epsw0 == 0.1 //done
																 //-1. <-????
}

double fa(double diffu, double A) {
	using namespace cont;
	return STIM11*(diffu > 0 ? diffu : 0) - A*Kaa;
}