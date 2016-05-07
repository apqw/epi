#pragma once
#include "define.h"
double _ljmain(double r1, double r2, double distSq);

double _ljmain_der_near(double r1, double r2, double dist);

double _ljmain_der_far(double r1, double r2, double dist);

double _adhe(double distlj, double rad_sum, double spring_const);

double min0(double u);

double fB(double age, double B, bool cornif_flg);

double fu(double u, double v, double p, double B);

double fw(double diff, double w);

double fa(double diffu, double A);