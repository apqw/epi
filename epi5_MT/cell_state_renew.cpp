#include "cell_state_renew.h"
#include "cell.h"
#include "cellmanager.h"
#include "define.h"
#include "Random_gen.h"
#include "utils.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <cassert>
inline void calc_dermis_normal(const Cell*const RESTRICT me, const Cell*const RESTRICT dermis, double& outx, double& outy, double& outz) {
	double nvx = p_diff_x(me->x(), dermis->x());
	double nvy = p_diff_y(me->y(), dermis->y());
	double nvz = me->z() - dermis->z();

	double norm = sqrt(DIST_SQ(nvx, nvy, nvz));
	outx = nvx / norm;
	outy = nvy / norm;
	outz = nvz / norm;

}

void div_direction(const Cell*const RESTRICT me, const Cell*const RESTRICT dermis, double* outx, double* outy, double* outz) {
	
	double nx, ny, nz;
	calc_dermis_normal(me, dermis, nx, ny, nz);
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

void divide_try(CellManager& cman, Cell*const RESTRICT div) {
	using namespace cont;
	if (div->pair != nullptr)return;
	const double div_gamma= stoch_corr_coef*DT_Cell*(div->is_malignant ? accel_div : 1)*eps_kb*ca2p_init / (div->div_age_thresh*stoch_div_time_ratio);

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
		div->ex_inert,
		0, 0,//set ages 0
		0, 0,//set fats 0
		delta_L,//init
		agki_max, //musume thresh
		div->rest_div_times,
		div->is_malignant
		);
	div->pair->pair = div;
	div->spr_nat_len=delta_L;
	div->ageb = 0;
	if (div->state == MUSUME) {
		div->rest_div_times--;
		div->pair->rest_div_times--;
	}

	assert(div->dermis() != nullptr);
	double divx, divy, divz;
	div_direction(div, div->dermis(), &divx, &divy, &divz);
	//!set value
	div->x._set(div->x()+ divx*0.5*delta_L);
	div->y._set(div->y()+divy*0.5*delta_L);
	div->z._set(div->z() + divz*0.5*delta_L);

	div->pair->x._set(div->pair->x()- divx*0.5*delta_L);
	div->pair->y._set(div->pair->y() - divy*0.5*delta_L);
	div->pair->z._set(div->pair->z() - divz*0.5*delta_L);
	printf("new cell detected.\n");
}

inline double ageb_const(const Cell*const RESTRICT c) {
	using namespace cont;
	return (c->is_malignant ? accel_div : 1)*eps_kb*(ca2p_init + alpha_b*min0(c->ca2p_avg - ca2p_init));
}

inline double agek_const(const Cell*const RESTRICT c) {
	using namespace cont;
	//assert(state() == DEAD || state() == AIR || state() == ALIVE);
	if (c->state==DEAD || c->state==AIR) {
		return eps_kk*ca2p_init;
	}
	else {
		return (c->is_malignant ? accel_diff : 1.0)*eps_ks*
			(S0 + alpha_k*min0(c->ca2p_avg - ca2p_init));
	}
}

inline double k_lipid_release(const Cell*const RESTRICT c) {
	using namespace cont;
	return 0.25*lipid_rel*(1 + tanh((c->ca2p_avg - ubar) / 0.01))*(1 + tanh((c->agek - THRESH_SP) / delta_lipid));
}

inline double k_lipid(const Cell*const RESTRICT c) {
	using namespace cont;
	return 0.25*lipid*(1 + tanh((ubar - c->ca2p_avg) / 0.01))*(1 + tanh((c->agek - THRESH_SP) / delta_lipid));
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

void _MUSUME_state_renew(CellManager& cman,Cell*const RESTRICT musume) {
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

	if (musume->rest_div_times > 0 && musume->ageb >= musume->div_age_thresh*(1.0 - cont::stoch_div_time_ratio)) {
		divide_try(cman, musume);
	}
	else {
		musume->ageb +=cont::DT_Cell*ageb_const(musume);
	}
}

void _FIX_state_renew(CellManager& cman, Cell*const RESTRICT fix) {
	if (fix->dermis() == nullptr) {
		printf("err\n");
		printf("x:%lf,y:%lf,z:%lf\n", fix->x(), fix->y(), fix->z());
        printf("connected_num:%zd\n", fix->connected_cell.size());
		assert(fix->dermis() != nullptr);
		return;
	}

	if (fix->ageb >= fix->div_age_thresh*(1.0 - cont::stoch_div_time_ratio)) {
		divide_try(cman, fix);
	}
	else {
		fix->ageb += cont::DT_Cell*ageb_const(fix);
	}
}

void _DEAD_AIR_state_renew(CellManager& cman, Cell*const RESTRICT da) {
	using namespace cont;
	if (da->agek >= ADHE_CONST&&da->connected_cell.size() <= DISA_conn_num_thresh) {
		da->state = DISA;
		cman.add_remove_queue(da->get_index());
	}
	else {
		da->agek += DT_Cell*agek_const(da);
	}
}

inline void _ALIVE_state_renew(CellManager& cman, Cell*const RESTRICT al) {
	using namespace cont;

	if (al->agek >= THRESH_DEAD) {
		cornificate(cman, al);
	}
	else {
		
		const double tmp = k_lipid_release(al)*al->in_fat;
		al->in_fat+=DT_Cell*(k_lipid(al)*(1.0 - al->in_fat) - tmp);
		al->ex_fat += DT_Cell*tmp;

		al->agek += DT_Cell*agek_const(al); //update last
	}
}

void pair_disperse(Cell*const RESTRICT c) {
	using namespace cont;
	assert(c->pair != nullptr);
	assert(c->pair->pair == c);
	double rad_sum = c->radius + c->pair->radius;
	double unpair_th = unpair_dist_coef*rad_sum;
	double distSq = 0;
	if (c->spr_nat_len < 2.0*c->radius) {
		c->spr_nat_len += DT_Cell*eps_L;
		c->pair->spr_nat_len = c->spr_nat_len;
	}
	else if ((distSq = p_cell_dist_sq(c, c->pair))>unpair_th*unpair_th) {
		c->spr_nat_len = 0;
		c->pair->spr_nat_len = 0;
		c->pair->pair = nullptr;
		c->pair = nullptr;
		printf("unpaired. distSq:%lf\n", distSq);
	}
}
///////////////////////////////////////////////////////////////////////

void cell_state_renew(CellManager & cman)
{
	cman.other_foreach_parallel_native([&](Cell*const RESTRICT c) {
		//auto&c = cman[i];
		switch (c->state) {
		case ALIVE:
			_ALIVE_state_renew(cman, c);
			break;
		case DEAD:case AIR:
			_DEAD_AIR_state_renew(cman, c);
			break;
		case FIX:
			_FIX_state_renew(cman, c);
			break;
		case MUSUME:
			_MUSUME_state_renew(cman, c);
			break;
		default:
			break;
		}
	});
	cman.remove_exec();
	cman.other_foreach_parallel_native([](Cell*const RESTRICT c) {
		//auto&c = cman[i];
		if (c->pair != nullptr&&no_double_count(c, c->pair)) {
			pair_disperse(c);
		}
	});
}

void initialize_sc(CellManager & cman,double zzmax)
{
	cman.other_foreach_parallel_native([zzmax](Cell*const RESTRICT c) {
		if (c->state == ALIVE&&zzmax - c->z() < 8 * cont::R_max) {
			c->agek = c->agek > cont::THRESH_DEAD ? c->agek : cont::THRESH_DEAD;
			c->state = AIR;
		}
	});
}
