#pragma once
#include "const.h"
#include "cell.h"
#include "field.h"

double fu(double u, double v, double p, double B)
{
	//  static double result_debug=0.0;

	return kf * (mu0 + mu1 * p / (p + kmu)) * (para_b + para_bb * u / (para_k1 + u)) * v
		- gamma * u / (kg + u) + beta + kbc * B*B * Cout / (Hb + B*B);

}
template<CELL_STATE state>
double th(double age)
{
	static_assert(state == ALIVE || state == FIX || state == MUSUME, "error in function fv");
	double tmp = 0;
	STATIC_IF(state == ALIVE) {
		tmp = thgra + ((thpri - thgra)*0.5) * (1.0 + tanh((THRESH_SP - age) / delta_th));
	}STATIC_ELSEIF(state == FIX || state == MUSUME) {
		tmp = thpri;
	}STATIC_ENDIF;
	return tmp;
}
template<CELL_STATE state>
double fv(const CellPtr& cell)
{
	double age = 0;
	STATIC_IF(state == MUSUME || state == FIX) {
		age = cell->old_data.ageb;
	}STATIC_ELSE{
		age = cell->old_data.agek;
	}STATIC_ENDIF;
	return (para_k2 * para_k2 / (para_k2 * para_k2 + cell->old_data.ca2p * cell->old_data.ca2p) - cell->old_data.ex_inert) / th<state>(age);
}

template<CELL_STATE state>
double Kpa(double age) {
	//static_assert(state == ALIVE || state == FIX || state == MUSUME, "error in function fv");
	double tmp = 0;
	STATIC_IF(state == ALIVE) {
		tmp = Kgra + ((Kpri - Kgra)*0.5) * (1.0 + tanh((THRESH_SP - age) / delta_K));
	}STATIC_ELSEIF(state == FIX || state == MUSUME) {
		tmp = Kpri;
	}STATIC_ENDIF;
	return tmp;
}

template<CELL_STATE state>
double h(double a, double age) {
	return Kpa<state>(age)*a / (H0 + a);
}

template<CELL_STATE state>
double fp(double a, const CellPtr& cell) {
	double age = 0;
	STATIC_IF(state == MUSUME || state == FIX) {
		age = cell->old_data.ageb;
	}STATIC_ELSE{
		age = cell->old_data.agek;
	}STATIC_ENDIF;
	return h<state>(a, age) - Kpp*cell->old_data.IP3;
}

template<CELL_STATE state>
double IAG(double age) {
	double tmp = 0;
	STATIC_IF(state == FIX || state == MUSUME) {
		tmp = iage_kitei;
	}STATIC_ELSE{
		tmp = 0.5*(1.0 + tanh((age - THRESH_SP) / delta_I));
	}STATIC_ENDIF;
	return tmp;
}

inline double fw(double diff, double w) {
	return (1. - w) + (-1. + tanh((wd - diff) / 0.1)) / 2.;
}

double fa(double diffu, double A) {
	return STIM11*(diffu > 0 ? diffu : 0) - A*Kaa;
}

double fB(double age, double B, int flg_cornified) {
	return flg_cornified*(age > THRESH_DEAD - DUR_ALIVE && age <= THRESH_DEAD + DUR_DEAD) - kb*B;
}

bool is_within_cell(double x1, double y1, double z1, double cx, double cy, double cz, double crad) {
	return dist3Sq(x1, y1, z1, cx, cy, cz) < crad*crad;
}

template<typename T>
inline T grid8_avg(const Arr3D<T>& grd, int ix, int iy, int iz) {
	return (grd[ix][iy][iz] + grd[ix + 1][iy][iz] + grd[ix][iy + 1][iz]
		+ grd[ix][iy][iz + 1] + grd[ix + 1][iy + 1][iz] + grd[ix + 1][iy][iz + 1]
		+ grd[ix][iy + 1][iz + 1] + grd[ix + 1][iy + 1][iz + 1]) / 8.0;
}