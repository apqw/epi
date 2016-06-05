#include "define.h"
#include "cell_interaction.h"
#include "cell.h"
#include "cellmanager.h"
#include "utils.h"
#include <functional>
#include <cmath>
#include <tbb/task_group.h>

#define CIFuncCName Coef
#define CIFuncDecl(name) struct name{static double CIFuncCName(const Cell* c1,const Cell* c2);}
template<class Fn>
void cell_interaction_apply(Cell* c1, Cell* c2) {
	double coef = Fn::CIFuncCName(c1, c2);
	double dumx = cont::DT_Cell*coef*p_diff_x(c1->x(), c2->x());
	double dumy = cont::DT_Cell*coef*p_diff_y(c1->y(), c2->y());
	double dumz = cont::DT_Cell*coef*(c1->z() - c2->z());
	c1->x += dumx;
	c1->y += dumy;
	c1->z += dumz;

	c2->x -= dumx;
	c2->y -= dumy;
	c2->z -= dumz;
}

CIFuncDecl(CI_der_to_der);
CIFuncDecl(CI_al_air_de_to_al_air_de_fix_mu);
CIFuncDecl(CI_fix_mu_to_fix_mu);
CIFuncDecl(CI_mu_to_memb);
CIFuncDecl(CI_fix_to_memb);
CIFuncDecl(CI_memb_to_memb);
CIFuncDecl(CI_other);
CIFuncDecl(CI_pair);

inline bool is_der_near(const Cell* der1, const Cell* der2) {
	return sqrt(p_cell_dist_sq(der1, der2))+ cont::delta_R < der1->radius + der2->radius;
}

inline double ljmain_der_near(const Cell* der1, const Cell* der2) {
	/*
		距離の2乗だけでは不可
	*/
	double dist = sqrt(p_cell_dist_sq(der1, der2));
	double dist2 = dist + cont::delta_R;
	double rad_sum = der1->radius + der2->radius;
	double LJ6 = rad_sum*rad_sum / (dist2*dist2);
	LJ6 = LJ6*LJ6*LJ6;
	return 4.0*cont::eps_m*LJ6*(LJ6 - 1) / (dist*dist2) + cont::para_ljp2;
}

inline bool is_near(const Cell* c1, const Cell* c2) {
	double rad_sum = c1->radius + c2->radius;
	return p_cell_dist_sq(c1, c2) < rad_sum*rad_sum;
}

inline double ljmain_der_far(const Cell* der1, const Cell* der2) {
	double dist = sqrt(p_cell_dist_sq(der1, der2));
	double dist2 = dist + cont::delta_R;
	double rad_sum = der1->radius + der2->radius;
	double LJ6 = rad_sum*rad_sum / (dist2*dist2);
	LJ6 = LJ6*LJ6*LJ6;
	return 4.0*cont::eps_m*LJ6*(LJ6 - 1) / (dist*dist) + cont::para_ljp2;
}

inline double ljmain(const Cell* c1, const Cell* c2) {
	double distSq = p_cell_dist_sq(c1, c2);
	double rad_sum = c1->radius + c2->radius;
	double LJ6 = rad_sum*rad_sum / distSq;
	LJ6 = LJ6*LJ6*LJ6;
	return 4.0*cont::eps_m*LJ6*(LJ6 - 1) / distSq;
}

inline double adhesion(const Cell* c1, const Cell* c2,double spring_const) {
	using namespace cont;
	double distlj = sqrt(p_cell_dist_sq(c1, c2));
	double rad_sum = c1->radius + c2->radius;
	double LJ2_m1;

	LJ2_m1 = (distlj / rad_sum) - 1;
	return (LJ2_m1 + 1 > LJ_THRESH ?
		0.0 :
		-(spring_const / distlj)
		* LJ2_m1*(1 - LJ2_m1*LJ2_m1 / ((LJ_THRESH - 1.0)*(LJ_THRESH - 1.0))));
}

/////////////////////////////////////////////////////////////////////////


double CI_der_to_der::CIFuncCName(const Cell* c1, const Cell* c2) {
	if (is_der_near(c1, c2)) {
		return ljmain_der_near(c1, c2);
	}
	else if (is_near(c1, c2)) {
		return cont::para_ljp2;
	}
	else {
		return ljmain_der_far(c1, c2);
	}
}

double CI_al_air_de_to_al_air_de_fix_mu::CIFuncCName(const Cell* c1, const Cell* c2) {
	if (is_near(c1, c2)) {
		return ljmain(c1, c2);
	}
	else {
		double spf = cont::K_DESMOSOME;
		if (c1->agek() > cont::THRESH_SP && c2->agek() > cont::THRESH_SP) {
			spf = cont::K_TOTAL;
		}
		return adhesion(c1, c2, spf);
	}
}

double CI_fix_mu_to_fix_mu::CIFuncCName(const Cell* c1, const Cell* c2) {
	//pairかどうかは外で判定しよう
	if (is_near(c1, c2)) {
		return ljmain(c1, c2);
	}
	else {
		return adhesion(c1, c2, cont::K_DESMOSOME);
	}
}

double CI_mu_to_memb::CIFuncCName(const Cell* c1, const Cell* c2) {
	//spring_force_to_membが設定されている必要がある
	if (is_near(c1, c2)) {
		return ljmain(c1, c2);
	}
	else {
		return adhesion(c1, c2, c1->spring_force_to_memb);
	}
}

double CI_fix_to_memb::CIFuncCName(const Cell* c1, const Cell* c2) {
	if (is_near(c1, c2)) {
		return ljmain(c1, c2);
	}
	else {
		double distlj = sqrt(p_cell_dist_sq(c1, c2));
		double LJ2 = distlj / (c1->radius + c2->radius);
		return -(cont::Kspring / distlj)*(LJ2 - 1.0);
	}
}

double CI_memb_to_memb::CIFuncCName(const Cell* c1, const Cell* c2) {
	using namespace cont;
	double rad_sum = c1->radius+ c2->radius;
	double cr_dist = rad_sum*P_MEMB;
	double lambda_dist = (1.0 + P_MEMB)*rad_sum;
	double distSq = p_cell_dist_sq(c1, c2);
	if (distSq< cr_dist*cr_dist) {
		//assert(fabs(ljmain(me, oppo)) < 1000);
		double LJ6 = cr_dist*cr_dist / distSq;
		LJ6 = LJ6*LJ6*LJ6;
		return 4.0*cont::eps_m*LJ6*(LJ6 - 1) / distSq;
	}
	else if (distSq < rad_sum*rad_sum) {
		double distlj = sqrt(distSq);
		return -(DER_DER_CONST / distlj) * (distlj / cr_dist - 1.0);
	}
	else {
		double distlj = sqrt(distSq);
		double LJ6 = cr_dist / (lambda_dist - distlj);
		LJ6 = LJ6*LJ6;
		LJ6 = LJ6*LJ6*LJ6;
		return -(DER_DER_CONST / rad_sum)*((1.0 - P_MEMB) / P_MEMB)
			- 4 * eps_m*(LJ6*(LJ6 - 1.0)) / ((lambda_dist - distlj)*distlj);
	}
}

double CI_other::CIFuncCName(const Cell* c1, const Cell* c2) {
	return ljmain(c1, c2);
}

double CI_pair::CIFuncCName(const Cell* c1, const Cell* c2) {
	double dist = sqrt(p_cell_dist_sq(c1, c2));
	double force = cont::Kspring_division*(dist - c1->spr_nat_len());
	return -force / dist;
}

void wall_interaction(Cell* c) {
	double distlj = 2.0*c->z();
	double LJ6 = c->radius / c->z();
	LJ6 = LJ6*LJ6;
	LJ6 = LJ6*LJ6*LJ6;
	double ljm = 4.0*cont::eps_m*LJ6*(LJ6 - 1.0) / (distlj*distlj);
	c->z += cont::DT_Cell* ljm*2.0*c->z();
}
////////////////////////////////////////////////////////////////////
void _memb_bend_calc1(Cell * memb)
{
	double dn = sqrt(p_cell_dist_sq(memb->md.memb_r, memb));

	memb->md.nv[0] = p_diff_x(memb->md.memb_r->x(), memb->x()) / dn;
	memb->md.nv[1] = p_diff_y(memb->md.memb_r->y(), memb->y()) / dn;
	memb->md.nv[2] = (memb->md.memb_r->z() - memb->z()) / dn;
	memb->md.dn = dn;

	double dm = sqrt(p_cell_dist_sq(memb->md.memb_u, memb));
	memb->md.mv[0] = p_diff_x(memb->md.memb_u->x(), memb->x()) / dm;
	memb->md.mv[1] = p_diff_y(memb->md.memb_u->y(), memb->y()) / dm;
	memb->md.mv[2] = (memb->md.memb_u->z() - memb->z()) / dm;
	memb->md.dm = dm;
}

void _memb_bend_calc2(Cell * memb)
{
	memb->md.ipn =
		memb->md.nv[0] * memb->md.memb_r->md.nv[0]
		+ memb->md.nv[1] * memb->md.memb_r->md.nv[1]
		+ memb->md.nv[2] * memb->md.memb_r->md.nv[2];

	memb->md.ipm =
		memb->md.mv[0] * memb->md.memb_u->md.mv[0]
		+ memb->md.mv[1] * memb->md.memb_u->md.mv[1]
		+ memb->md.mv[2] * memb->md.memb_u->md.mv[2];
}

void _memb_bend_calc3(Cell * memb)
{
	auto& nvr = memb->md.memb_r->md.nv;
	auto& nvl = memb->md.memb_l->md.nv;
	auto& nvll = memb->md.memb_ll->md.nv;

	auto& mvu = memb->md.memb_u->md.mv;
	auto& mvb = memb->md.memb_b->md.mv;
	auto& mvbb = memb->md.memb_bb->md.mv;

	auto& ipnl = memb->md.memb_l->md.ipn;
	auto& ipnll = memb->md.memb_ll->md.ipn;

	auto& ipmb = memb->md.memb_b->md.ipm;
	auto& ipmbb = memb->md.memb_bb->md.ipn;

	auto& dnl = memb->md.memb_l->md.dn;
	auto& dmb = memb->md.memb_b->md.dm;

	auto& dn = memb->md.dn;
	auto& dm = memb->md.dm;
	auto& nv = memb->md.nv;
	auto& mv = memb->md.mv;
	auto& ipn = memb->md.ipn;
	auto& ipm = memb->md.ipm;

	double x = cont::KBEND*(
		-(1.0 - ipn)*(nvr[0] - ipn*nv[0]) / dn
		+ (1.0 - ipnl)*((nv[0] - ipnl*nvl[0]) / dnl - (nvl[0] - ipnl*nv[0]) / dn)
		+ (1.0 - ipnll)*(nvll[0] - ipnll*nvl[0]) / dnl

		- (1.0 - ipm)*(mvu[0] - ipm*mv[0]) / dm
		+ (1.0 - ipmb)*((mv[0] - ipmb*mvb[0]) / dmb - (mvb[0] - ipmb*mv[0]) / dm)
		+ (1.0 - ipmbb)*(mvbb[0] - ipmbb*mvb[0]) / dmb);

	double y = cont::KBEND*(
		-(1.0 - ipn)*(nvr[1] - ipn*nv[1]) / dn
		+ (1.0 - ipnl)*((nv[1] - ipnl*nvl[1]) / dnl - (nvl[1] - ipnl*nv[1]) / dn)
		+ (1.0 - ipnll)*(nvll[1] - ipnll*nvl[1]) / dnl

		- (1.0 - ipm)*(mvu[1] - ipm*mv[1]) / dm
		+ (1.0 - ipmb)*((mv[1] - ipmb*mvb[1]) / dmb - (mvb[1] - ipmb*mv[1]) / dm)
		+ (1.0 - ipmbb)*(mvbb[1] - ipmbb*mvb[1]) / dmb);

	double z = cont::KBEND*(
		-(1.0 - ipn)*(nvr[2] - ipn*nv[2]) / dn
		+ (1.0 - ipnl)*((nv[2] - ipnl*nvl[2]) / dnl - (nvl[2] - ipnl*nv[2]) / dn)
		+ (1.0 - ipnll)*(nvll[2] - ipnll*nvl[2]) / dnl

		- (1.0 - ipm)*(mvu[2] - ipm*mv[2]) / dm
		+ (1.0 - ipmb)*((mv[2] - ipmb*mvb[2]) / dmb - (mvb[2] - ipmb*mv[2]) / dm)
		+ (1.0 - ipmbb)*(mvbb[2] - ipmbb*mvb[2]) / dmb);

	memb->x += cont::DT_Cell*x;
	memb->y += cont::DT_Cell*y;
	memb->z += cont::DT_Cell*z;
}

void memb_bend(CellManager& cman) {
	cman.memb_foreach_parallel([&](size_t i) {
		_memb_bend_calc1(cman[i]);
	});

	cman.memb_foreach_parallel([&](size_t i) {
		_memb_bend_calc2(cman[i]);
	});

	cman.memb_foreach_parallel([&](size_t i) {
		_memb_bend_calc3(cman[i]);
	});
}

//////////////////////////////////////////////////////////////////////////////////////

inline bool collide_with_wall(Cell* c) {
	return c->z() < c->radius;
}



void _MEMB_interaction(Cell* memb) {
	if (collide_with_wall(memb)) {
		wall_interaction(memb);
	}
	memb->connected_cell.foreach([&](Cell* conn) {
		if (conn->state == MEMB && no_double_count(memb,conn)) {
			cell_interaction_apply<CI_memb_to_memb>(memb, conn);
		}
	});
}

void _DER_interaction(Cell* der) {
	if (collide_with_wall(der)) {
		wall_interaction(der);
	}
	der->connected_cell.foreach([&](Cell* conn) {
		if (conn->state == DER && no_double_count(der, conn)) {
			cell_interaction_apply<CI_der_to_der>(der, conn);
		}
		else {
			cell_interaction_apply<CI_other>(der, conn);
		}
	});
}

void _AL_AIR_DE_interaction(Cell* aad) {
	aad->connected_cell.foreach([&](Cell* conn) {
		switch (conn->state) {
		case ALIVE:case AIR:case DEAD:
			if (no_double_count(aad, conn)) {
				cell_interaction_apply<CI_al_air_de_to_al_air_de_fix_mu>(aad, conn);
			}
			break;
		case DER:
			break;
		default:
			cell_interaction_apply<CI_al_air_de_to_al_air_de_fix_mu>(aad, conn);
			break;
		}
	
	});
}

void _FIX_interaction(Cell* fix) {
	fix->connected_cell.foreach([&](Cell* conn) {
		switch (conn->state) {
		case FIX:
			if (fix->pair!=conn&&no_double_count(fix, conn)) {
				cell_interaction_apply<CI_fix_mu_to_fix_mu>(fix, conn);
			}
			break;
		case MUSUME:
			if(fix->pair!=conn)cell_interaction_apply<CI_fix_mu_to_fix_mu>(fix, conn);
			break;
		case MEMB:
			if (fix->dermis() == conn) {
				cell_interaction_apply<CI_fix_to_memb>(fix, conn);
			}
			else {
				cell_interaction_apply<CI_other>(fix, conn);
			}
			break;
		default:
			break;
		}
	});
}


bool paired_with_fix(Cell* c) {
	return c->pair != nullptr&&c->pair->state == FIX;
}
void _MUSUME_interaction(Cell* musume) {
	musume->spring_force_to_memb = 0;
	if (paired_with_fix(musume)) {
		musume->spring_force_to_memb = cont::Kspring;
	}
	else if (musume->rest_div_times > 0) {
		musume->spring_force_to_memb = cont::Kspring_d;
	}
	musume->connected_cell.foreach([&](Cell* conn) {
		if (conn->state == MEMB) {
			if (musume->dermis() == conn) {
				cell_interaction_apply<CI_mu_to_memb>(musume, conn);
			}
			else {
				cell_interaction_apply<CI_other>(musume, conn);
			}
		}
		else if (conn->state == MUSUME &&musume->pair!=conn&& no_double_count(musume, conn)) {
			cell_interaction_apply<CI_fix_mu_to_fix_mu>(musume, conn);
		}
	});
}

void _NULL_interaction(Cell* n) {

}

void _pair_interaction(Cell* paired) {
	if (no_double_count(paired, paired->pair)) {
		cell_interaction_apply<CI_pair>(paired, paired->pair);
	}
}

void(*interaction_table[cont::STATE_NUM])(Cell*);

void _init_interaction_table(void(*i_tbl[cont::STATE_NUM])(Cell*)) {
	i_tbl[ALIVE] = _AL_AIR_DE_interaction;
	i_tbl[DEAD] = _AL_AIR_DE_interaction;
	i_tbl[DISA] = _NULL_interaction;
	i_tbl[UNUSED] = _NULL_interaction;
	i_tbl[FIX] = _FIX_interaction;
	i_tbl[BLANK] = _NULL_interaction;
	i_tbl[DER] = _DER_interaction;
	i_tbl[MUSUME] = _MUSUME_interaction;
	i_tbl[AIR] = _AL_AIR_DE_interaction;
	i_tbl[MEMB] = _MEMB_interaction;
}


//////////////////////////////////////////////////////////////////////////////////////


void init_cell_interaction() {
	_init_interaction_table(interaction_table);
}

void cell_interaction(CellManager & cman)
{
	tbb::task_group t;
	t.run([&] {
		memb_bend(cman);
	});
	t.run([&] {
		cman.all_foreach_parallel([&](size_t i) {
			auto&c = cman[i];
			interaction_table[c->state](c);
			if (c->pair != nullptr)_pair_interaction(c);
		});
	});
	t.wait();
		

	
}
