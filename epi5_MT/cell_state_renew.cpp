#include "cell_state_renew.h"
#include "cell.h"
#include "cellmanager.h"
#include "define.h"
#include "Random_gen.h"
#include "utils.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <cassert>
void calc_dermis_normal(const Cell* me, const Cell*dermis, double* outx, double* outy, double* outz) {
	double nvx = p_diff_x(me->x(), dermis->x());
	double nvy = p_diff_y(me->y(), dermis->y());
	double nvz = me->z() - dermis->z();

	double norm = sqrt(DIST_SQ(nvx, nvy, nvz));
	*outx = nvx / norm;
	*outy = nvy / norm;
	*outz = nvz / norm;

}

void div_direction(const Cell* me, const Cell*dermis, double* outx, double* outy, double* outz) {
	
	double nx, ny, nz;
	calc_dermis_normal(me, dermis, &nx, &ny, &nz);
	double ox, oy, oz;
	double sum;
	do {
		double rand_theta = M_PI*genrand_real();
		double cr1 = cos(rand_theta);
		double sr1 = sin(rand_theta);
		double rand_phi = 2 * M_PI*genrand_real();
		double cr2 = cos(rand_phi);
		double sr2 = sin(rand_phi);
		double rvx = sr1*cr2;
		double rvy = sr1*sr2;
		double rvz = cr1;

		ox = ny*rvz - nz*rvy;
		oy = nz*rvx - nx*rvz;
		oz = nx*rvy - ny*rvx;
	} while ((sum = DIST_SQ(ox, oy, oz)) < 1.0e-14);
	sum = sqrt(sum);
	*outx = ox / sum;
	*outy = oy / sum;
	*outz = oz / sum;
}

void divide_try(CellManager& cman, Cell* div) {
	using namespace cont;
	if (div->pair != nullptr)return;
	double div_gamma= stoch_corr_coef*DT_Cell*(div->is_malignant ? accel_div : 1)*eps_kb*ca2p_init / (div->div_age_thresh*stoch_div_time_ratio);

	if (STOCHASTIC&&genrand_real() > div_gamma) {
		return ;
	}

	div->pair = cman.create(
		MUSUME,
		div->x(), div->y(), div->z(),
		div->radius,
		div->ca2p(),
		div->ca2p(),//this is avg value,do not use orignal avg
		div->IP3(),
		div->ex_inert(),
		0, 0,//set ages 0
		0, 0,//set fats 0
		delta_L,//init
		agki_max, //musume thresh
		div->rest_div_times,
		div->is_malignant
		);
	div->pair->pair = div;
	div->spr_nat_len._set(delta_L);
	div->ageb._set(0);
	if (div->state == MUSUME) {
		div->rest_div_times--;
		div->pair->rest_div_times--;
	}

	assert(div->dermis() != nullptr);
	double divx, divy, divz;
	div_direction(div, div->dermis(), &divx, &divy, &divz);
	div->x += divx*0.5*delta_L;
	div->y += divy*0.5*delta_L;
	div->z += divz*0.5*delta_L;

	div->pair->x -= divx*0.5*delta_L;
	div->pair->y -= divy*0.5*delta_L;
	div->pair->z -= divz*0.5*delta_L;
}

double ageb_const(Cell* c) {
	using namespace cont;
	return (c->is_malignant ? accel_div : 1)*eps_kb*(ca2p_init + alpha_b*min0(c->ca2p_avg - ca2p_init));
}

double agek_const(Cell* c) {
	using namespace cont;
	//assert(state() == DEAD || state() == AIR || state() == ALIVE);
	if (c->state == DEAD || c->state == AIR) {
		return eps_kk*ca2p_init;
	}
	else {
		return (c->is_malignant ? accel_diff : 1.0)*eps_ks*
			(S0 + alpha_k*min0(c->ca2p_avg - ca2p_init));
	}
}

double k_lipid_release(Cell* c) {
	using namespace cont;
	return 0.25*lipid_rel*(1 + tanh((c->ca2p_avg - ubar) / 0.01))*(1 + tanh((c->agek() - THRESH_SP) / delta_lipid));
}

double k_lipid(Cell* c) {
	using namespace cont;
	return 0.25*lipid*(1 + tanh((ubar - c->ca2p_avg) / 0.01))*(1 + tanh((c->agek() - THRESH_SP) / delta_lipid));
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

void _MUSUME_state_renew(CellManager& cman,Cell* musume) {
	if (musume->dermis() == nullptr && musume->pair == nullptr) {
		if (SYSTEM == WHOLE) {
			printf("ALIVE detected\n");
			musume->state = ALIVE;
		}else if (SYSTEM == BASAL) {
			musume->state = DISA;
			cman.add_remove_queue(musume->get_index());
		}
		return;
	}

	if (musume->rest_div_times > 0 && musume->ageb() >= musume->div_age_thresh*(1.0 - cont::stoch_div_time_ratio)) {
		divide_try(cman, musume);
	}
	else {
		musume->ageb = musume->ageb()+cont::DT_Cell*ageb_const(musume);
	}
}

void _FIX_state_renew(CellManager& cman, Cell* fix) {
	if (fix->dermis() == nullptr) {
		printf("err\n");
		printf("x:%lf,y:%lf,z:%lf\n", fix->x(), fix->y(), fix->z());
		printf("connected_num:%d\n", fix->connected_cell.size());
		assert(fix->dermis() != nullptr);
		return;
	}

	if (fix->ageb() >= fix->div_age_thresh*(1.0 - cont::stoch_div_time_ratio)) {
		divide_try(cman, fix);
	}
	else {
		fix->ageb =fix->ageb()+ cont::DT_Cell*ageb_const(fix);
	}
}

void _DEAD_AIR_state_renew(CellManager& cman, Cell* da) {
	using namespace cont;
	if (da->agek() >= ADHE_CONST&&da->connected_cell.size() <= DISA_conn_num_thresh) {
		da->state = DISA;
		cman.add_remove_queue(da->get_index());
	}
	else {
		da->agek =da->agek()+ DT_Cell*agek_const(da);
	}
}

void _ALIVE_state_renew(CellManager& cman, Cell* al) {
	using namespace cont;

	if (al->agek() >= THRESH_DEAD) {
		al->state = DEAD;
		printf("sw updated:%d\n", ++cman.sw);
	}
	else {
		al->agek =al->agek()+ DT_Cell*agek_const(al);
		double tmp = k_lipid_release(al)*al->in_fat();
		al->in_fat=al->in_fat()+DT_Cell*(k_lipid(al)*(1.0 - al->in_fat()) - tmp);
		al->ex_fat =al->ex_fat()+ DT_Cell*tmp;
	}
}

void _NULL_state_renew(CellManager& cman, Cell* al) {

}

void(*state_renew_table[cont::STATE_NUM])(CellManager&,Cell*);

void _init_interaction_table(void(*sr_tbl[cont::STATE_NUM])(CellManager&, Cell*)) {
	sr_tbl[ALIVE] = _ALIVE_state_renew;
	sr_tbl[DEAD] = _DEAD_AIR_state_renew;
	sr_tbl[DISA] = _NULL_state_renew;
	sr_tbl[UNUSED] = _NULL_state_renew;
	sr_tbl[FIX] = _FIX_state_renew;
	sr_tbl[BLANK] = _NULL_state_renew;
	sr_tbl[DER] = _NULL_state_renew;
	sr_tbl[MUSUME] = _MUSUME_state_renew;
	sr_tbl[AIR] = _DEAD_AIR_state_renew;
	sr_tbl[MEMB] = _NULL_state_renew;
}


void pair_disperse(Cell* c) {
	using namespace cont;
	assert(c->pair != nullptr);
	assert(c->pair->pair == c);
	double rad_sum = c->radius + c->pair->radius;
	double unpair_th = unpair_dist_coef*rad_sum;
	double distSq = 0;
	if (c->spr_nat_len() < 2.0*c->radius) {
		c->spr_nat_len =c->spr_nat_len()+ DT_Cell*eps_L;
		c->pair->spr_nat_len._set(c->spr_nat_len() + DT_Cell*eps_L);
	}
	else if ((distSq = p_cell_dist_sq(c, c->pair))>unpair_th*unpair_th) {
		c->spr_nat_len._set(0);
		c->pair->spr_nat_len._set(0);
		c->pair->pair = nullptr;
		c->pair = nullptr;
		printf("unpaired. distSq:%lf\n", distSq);
	}
}
///////////////////////////////////////////////////////////////////////

void init_cell_state_renewal() {
	_init_interaction_table(state_renew_table);
}

void cell_state_renew(CellManager & cman)
{
	cman.other_foreach_parallel([&](size_t i) {
		auto&c = cman[i];
		state_renew_table[c->state](cman, c);
	});
	cman.remove_exec();
	cman.other_foreach_parallel([&](size_t i) {
		auto&c = cman[i];
		if (c->pair != nullptr&&no_double_count(c, c->pair)) {
			pair_disperse(c);
		}
	});
}
