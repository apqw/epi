#pragma once
#include <memory>
#include <cassert>
#include <unordered_map>


#include "const.h"
#include "funcs.h"
#include "static_if.hpp"

class Cell; //prototype decl
using CellPtr = std::shared_ptr<Cell>;
using CellIFunc = void(*)(CellPtr&, CellPtr&);
using CellICoefFunc = double(*)(const CellPtr&, const CellPtr&, double, double, double);
using CellID = uint32_t;

struct CellData {
	double x; double y; double z;

	double radius;

	double ca2p, ca2p_avg;

	double IP3;

	double agek, ageb;

	double re_ex_fat, in_fat;

	double pair_spr_nat_len;

	double spring_force;

	double div_age_thresh;

	double poisson_div_thresh;

	int rest_div_times;

	bool infinite_div;
	bool malignant;
	bool touch;

	const CellID src_stem_id;

	CellPtr div_pair_cell;
	CellPtr dermis;

	Arr2D<CellPtr> connected_cell;
	Arr2D<double> gj;
	double ex_inert;
	int connected_count;

	//std::unordered_map<CellID, double> ca2p_inert;

	CellData(
		uint32_t _src_stem_id = 0,
		double _x = 0,
		double _y = 0,
		double _z = 0,
		double _radius = 0,
		double _ca2p = u0,
		double _ca2p_avg = u0,
		double _re_ex_fat = 0,
		double _in_fat = 0,
		double _pair_spr_nat_len = 0,
		int _rest_div_times = 0,
		double _agek = 0,
		double _ageb = 0,
		double _IP3 = p0,
		double _spring_force = 0,
		double _div_age = 0,
		double _ex_inert = v0
		) :src_stem_id(_src_stem_id), x(_x), y(_y), z(_z), radius(_radius), ca2p(_ca2p), ca2p_avg(_ca2p_avg), re_ex_fat(_re_ex_fat), in_fat(_in_fat),
		pair_spr_nat_len(_pair_spr_nat_len), rest_div_times(_rest_div_times),
		agek(_agek), ageb(_ageb), IP3(_IP3), div_pair_cell(nullptr), dermis(nullptr), spring_force(_spring_force),
		div_age_thresh(_div_age), infinite_div(false), malignant(false),ex_inert(_ex_inert) {
		gj = Arr2D<double>(STATE_NUM, Arr1D<double>(0));
		connected_cell = Arr2D<CellPtr>(STATE_NUM, Arr1D<CellPtr>(0));
		for (auto& cc : connected_cell) {
			cc.reserve(N2);
		}

	};
	CellData& operator=(const CellData&ot) {

		assert(src_stem_id == ot.src_stem_id);

		x = ot.x;
		y = ot.y;
		z = ot.z;
		radius = ot.radius;
		ca2p = ot.ca2p;
		ca2p_avg = ot.ca2p_avg;
		re_ex_fat = ot.re_ex_fat;
		in_fat = ot.in_fat;
		pair_spr_nat_len = ot.pair_spr_nat_len;
		rest_div_times = ot.rest_div_times;
		agek = ot.agek;
		ageb = ot.ageb;
		IP3 = ot.IP3;
		ex_inert = ot.ex_inert;
		connected_cell = ot.connected_cell;//copied
		//ca2p_inert = ot.ca2p_inert;
		//div_pair_cell = ot.div_pair_cell;
		return *this;
	}
};

class Cell {
private:
	void div_direction(const CellData& me_data, const CellData& dermis_data, double *nx, double *ny, double *nz);
public:
	CellData data;
	CellData old_data;
	CellData div_cell_data;
	const CellID id;
	bool pending_kill = false;
	bool pending_born = true;
	bool state_changed = false;
	bool got_unpaired = false;
	CELL_STATE state_dest = UNUSED;

	

	void copy_data_to_old() {
		old_data = data;
	}

	//need dermis on new position
	template<CELL_STATE fix_or_musume>
	bool check_divide(CellPtr& my_shared_ptr) {
		//位置は更新済みなので注意
		static_assert(fix_or_musume == FIX || fix_or_musume == MUSUME, "non-fix_or_musume cant divide");
		if (old_data.rest_div_times <= 0
			|| old_data.ageb < old_data.div_age_thresh*(1 - stoch_div_time_ratio))return false;
		if (old_data.div_pair_cell != nullptr)return false;
		double div_th = DT_Cell*eps_kb*u0 / (old_data.div_age_thresh*stoch_div_time_ratio);
		if (old_data.malignant)div_th *= accel_div;

		if (STOCHASTIC && genrand_real() > div_th)return false;
		//divide!
		//set div cell data
		double orig_x = old_data.x;
		double orig_y = old_data.y;
		double orig_z = old_data.z;

		div_cell_data = old_data; //位置だけ更新済み

		div_cell_data.div_pair_cell = my_shared_ptr; //important

		if (!old_data.infinite_div) {
			div_cell_data.rest_div_times--;
			data.rest_div_times = old_data.rest_div_times - 1;
		}
		double nx, ny, nz;
		assert(old_data.dermis != nullptr);
		div_direction(old_data, old_data.dermis->old_data, &nx, &ny, &nz);

		double rh = 0.5*delta_L;
		data.x = orig_x + nx*rh;
		data.y = orig_y + ny*rh;
		data.z = orig_z + nz*rh;

		div_cell_data.x = orig_x - nx*rh;
		div_cell_data.y = orig_y - ny*rh;
		div_cell_data.z = orig_z - nz*rh;


		data.ageb = div_cell_data.ageb = 0;
		data.spring_force = div_cell_data.spring_force = 0;
		STATIC_IF(fix_or_musume == FIX) {
			data.spring_force = div_cell_data.spring_force = Kspring;
		}STATIC_ELSEIF(fix_or_musume == MUSUME) {
			data.spring_force = div_cell_data.spring_force = Kspring_d;
		}STATIC_ENDIF;
		div_cell_data.agek = 0;
		div_cell_data.re_ex_fat = 0;
		div_cell_data.in_fat = 0;

		return true;

	}

	void set_change_state(CELL_STATE cst) {
		assert(!state_changed); //no double change
		state_dest = cst;
		state_changed = true;
	}

	void mark_as_delete() {
		assert(!pending_kill); //no double kill
		pending_kill = true;
	}

	void mark_as_got_unpaired() {
		assert(!got_unpaired); //no double kill
		got_unpaired = true;
	}

	Cell(const CellID _id) :id(_id), pending_kill(false), pending_born(true), state_changed(false)
		, state_dest(UNUSED), got_unpaired(false) {};

	Cell(const CellID _id, const CellData &_data) :id(_id), data(_data), old_data(_data), pending_kill(false)
		, pending_born(true), state_changed(false), state_dest(UNUSED), got_unpaired(false) {};
};
template<class f>
void cell_intr_tmpl(CellPtr& me, CellPtr& oppo) {
	double mex = me->old_data.x;
	double mey = me->old_data.y;
	double mez = me->old_data.z;

	double ox = oppo->old_data.x;
	double oy = oppo->old_data.y;
	double oz = oppo->old_data.z;

	double diffx = perio_diff_x(mex, ox);
	double diffy = perio_diff_y(mey, oy);
	double diffz = mez - oz;

	double dum_ljmain = f()(me, oppo, diffx, diffy, diffz); //optimization may work?

	double dumxx = DT_Cell * dum_ljmain * diffx;
	double dumyy = DT_Cell * dum_ljmain * diffy;
	double dumzz = DT_Cell * dum_ljmain * diffz;
	assert(!isnan(dumxx));
	me->data.x = mex + dumxx;
	me->data.y = mey + dumyy;
	me->data.z = mez + dumzz;

	oppo->data.x = ox - dumxx;
	oppo->data.y = oy - dumyy;
	oppo->data.z = oz - dumzz;
}

void interac_wall(CellPtr& me);
class memb_memb_coef:FObj {
public:
	double operator()(const CellPtr& me, const CellPtr& oppo, double diffx, double diffy, double diffz);
};
class null_coef :FObj {
public:
	double operator()(const CellPtr& me, const CellPtr& oppo, double diffx, double diffy, double diffz);
};
class common_common_coef :FObj {
public:
	double operator()(const CellPtr& me, const CellPtr& oppo, double diffx, double diffy, double diffz);
};
class der_der_coef :FObj {
public:
	double operator()(const CellPtr& me, const CellPtr& oppo, double diffx, double diffy, double diffz);
};

class supra_others_coef :FObj {
public:
	double operator()(const CellPtr& me, const CellPtr& oppo, double diffx, double diffy, double diffz);
};
template<CELL_STATE fix_or_musume>
class stem_memb_coef :FObj {
public:
	double operator()(const CellPtr& me, const CellPtr& oppo, double diffx, double diffy, double diffz) {
		static_assert(fix_or_musume == FIX || fix_or_musume == MUSUME, "non-stem cell specified in stem_memb_coef");
		if (me->old_data.dermis == oppo) {

			double rad_sum = me->old_data.radius + oppo->old_data.radius;
			double rad_sumSq = rad_sum*rad_sum;
			double distSq = normSq(diffx, diffy, diffz);
			if (distSq < rad_sumSq) {
				double LJ1Sq = rad_sumSq / distSq;
				double LJ6 = LJ1Sq*LJ1Sq*LJ1Sq;
				return 4.0*eps_m * LJ6 * (LJ6 - 1.0) / distSq;
			}
			else {

				double tmp = 0;
				STATIC_IF(fix_or_musume == MUSUME) {
					tmp = adhesion(sqrt(distSq), rad_sum, me->old_data.spring_force);
				}STATIC_ELSE{
					double dist = sqrt(distSq);
				tmp = -(me->old_data.spring_force / dist)*(dist / rad_sum - 1.0);
				}STATIC_ENDIF;
				return tmp;
			}
		}
		else {
			return common_common_coef()(me, oppo, diffx, diffy, diffz);
		}
	}
};

class stem_stem_coef :FObj {
public:
	double operator()(const CellPtr& me, const CellPtr& oppo, double diffx, double diffy, double diffz);
};
class pair_coef :FObj {
public:
	double operator()(const CellPtr& me, const CellPtr& oppo, double diffx, double diffy, double diffz);
};
template<CELL_STATE state>
double ageb_const(const CellPtr& me) {
	static_assert(state == FIX || state == MUSUME, "ageb_const is for FIX and MUSUME");
	double cadif = me->old_data.ca2p_avg - u0;

	double tmp = eps_kb*(u0 + alpha_b*(cadif > 0 ? cadif : 0));
	if (me->old_data.malignant) {
		tmp *= accel_div;
	}
	return tmp;
}

template<CELL_STATE state>
double agek_const(const CellPtr& me) {
	static_assert(state == DEAD || state == AIR || state == ALIVE, "error: state must be DEAD (AIR) or ALIVE\n");
	double ret = 0;
	STATIC_IF(state == DEAD || state == AIR) {
		ret = eps_kk*u0;
	}STATIC_ELSEIF(state == ALIVE) {
		double diffu = me->old_data.ca2p_avg - u0;
		ret = eps_ks*(S0 + alpha_k*(diffu > 0 ? diffu : 0));
		if (me->old_data.malignant) {
			ret *= accel_diff;
		}
	}STATIC_ENDIF;
	return ret;
}

double k_lipid_release(double u, double age);

double k_lipid(double u, double age);

CellPtr find_dermis(const CellPtr& me);

template<bool with_new = false>
inline double cell_distSq(const CellPtr& c1, const CellPtr& c2) {
	double tmp;
	STATIC_IF(with_new) {
		tmp = dist3Sq(c1->data.x, c1->data.y, c1->data.z,
			c2->data.x, c2->data.y, c2->data.z);
	}STATIC_ELSE{
		tmp = dist3Sq(c1->old_data.x, c1->old_data.y, c1->old_data.z,
		c2->old_data.x, c2->old_data.y, c2->old_data.z);
	}STATIC_ENDIF;
	return tmp;
}
