#include "define.h"
#include "cell_interaction.h"
#include "cell.h"
#include "cellmanager.h"
#include "utils.h"
#include <functional>
#include <cmath>
#include <tbb/task_group.h>
/**
 *  @file ×–EŠÔ‚Ì‘ŠŒÝì—p(Žå‚ÉˆÊ’u)‚ÉŠÖ‚·‚é’è‹`
 */

#define POW6(x) ((x)*(x)*(x)*(x)*(x)*(x))
#define POW3(x) ((x)*(x)*(x))
#define CIFuncCName Coef
#define CIFuncDecl(name) struct name{static double CIFuncCName(const Cell*const RESTRICT c1,const Cell*const RESTRICT c2);}

/** LJƒ|ƒeƒ“ƒVƒƒƒ‹‚ÌŒW” */
static constexpr double eps_m				= 0.01;

/*
	MEMB_calc
*/
static constexpr double P_MEMB				= 1. / cont::COMPRESS_FACTOR;

/** L‚Ñ’e«ŒW” */
static constexpr double DER_DER_CONST		= 0.05;
static constexpr double K_TOTAL				= 3.0;
static constexpr double K_DESMOSOME_RATIO	= 0.01;
static constexpr double K_DESMOSOME			= K_TOTAL*K_DESMOSOME_RATIO;
static constexpr double Kspring				= 25.0;
static constexpr double Kspring_d			= 5.0;

/** ‹È‚°’e«ŒW” */
static constexpr double KBEND = 0.25;//0.5->0.25

/*
	DER_calc
*/
static constexpr double delta_R				= 0.4*cont::R_der;

static constexpr double para_ljp2			= 0.005;
static constexpr double Kspring_division	= 5.0;


template<class Fn>
inline void cell_interaction_apply(Cell*const RESTRICT c1, Cell*const RESTRICT c2) {
	const double coef = Fn::CIFuncCName(c1, c2);
    if (fabs(coef) > 100) {
        printf("Too strong interaction:%d and %d\n", (int)(c1->state), (int)(c2->state));
        assert(false);
        exit(1);
    }
	const double dumx = cont::DT_Cell*coef*p_diff_x(c1->x(), c2->x());
	const double dumy = cont::DT_Cell*coef*p_diff_y(c1->y(), c2->y());
	const double dumz = cont::DT_Cell*coef*(c1->z() - c2->z());
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
CIFuncDecl(CI_memb_to_memb_diag);
inline bool is_der_near(const Cell*const RESTRICT der1, const Cell*const RESTRICT der2) {
	return sqrt(p_cell_dist_sq(der1, der2))+ delta_R < der1->radius + der2->radius;
}

inline double ljmain_der_near(const Cell*const RESTRICT der1, const Cell*const RESTRICT der2) {
	/*
		‹——£‚Ì2æ‚¾‚¯‚Å‚Í•s‰Â
	*/
	const double dist = sqrt(p_cell_dist_sq(der1, der2));
	const double dist2 = dist + delta_R;
	const double rad_sum = der1->radius + der2->radius;
	const double LJ6 = POW6(rad_sum) / POW6(dist2);
	//LJ6 = LJ6*LJ6*LJ6;
	return 4.0*eps_m*LJ6*(LJ6 - 1) / (dist*dist2) + para_ljp2;
}

inline bool is_near(const Cell*const RESTRICT c1, const Cell*const RESTRICT c2) {
	//double rad_sum = c1->radius + c2->radius;
	return p_cell_dist_sq(c1, c2) < (c1->radius + c2->radius)*(c1->radius + c2->radius);
}

inline double ljmain_der_far(const Cell*const RESTRICT der1, const Cell*const RESTRICT der2) {
	const double dist = sqrt(p_cell_dist_sq(der1, der2));
	const double dist2 = dist + delta_R;
	const double rad_sum = der1->radius + der2->radius;
	const double LJ6 = POW6(rad_sum) / POW6(dist2);
	return 4.0*eps_m*LJ6*(LJ6 - 1) / (dist*dist) + para_ljp2;
}

inline double ljmain(const Cell*const RESTRICT c1, const  Cell*const RESTRICT c2) {

	const double distSq = p_cell_dist_sq(c1, c2);
	const double rad_sum = c1->radius + c2->radius;
	const double LJ6 = POW6(rad_sum) / POW3(distSq);
	//LJ6 = LJ6*LJ6*LJ6;
	return 4.0*eps_m*LJ6*(LJ6 - 1) / distSq;

}

inline double adhesion(const Cell*const RESTRICT c1, const  Cell*const RESTRICT c2,const double spring_const) {
	using namespace cont;
	const double distlj = sqrt(p_cell_dist_sq(c1, c2));
	const double rad_sum = c1->radius + c2->radius;
	//double LJ2_m1;

	const double LJ2_m1 = (distlj / rad_sum) - 1;
	return (LJ2_m1 + 1 > LJ_THRESH ?
		0.0 :
		-(spring_const / distlj)
		* LJ2_m1*(1 - LJ2_m1*LJ2_m1 / ((LJ_THRESH - 1.0)*(LJ_THRESH - 1.0))));
}

/////////////////////////////////////////////////////////////////////////


double CI_der_to_der::CIFuncCName(const Cell*const RESTRICT c1, const Cell*const RESTRICT c2) {
	if (is_der_near(c1, c2)) {
		return ljmain_der_near(c1, c2);
	}
	else if (is_near(c1, c2)) {
		return para_ljp2;
	}
	else {
		return ljmain_der_far(c1, c2);
	}
}

double CI_al_air_de_to_al_air_de_fix_mu::CIFuncCName(const Cell*const RESTRICT c1, const Cell*const RESTRICT c2) {
	if (is_near(c1, c2)) {
		return ljmain(c1, c2);
	}
	else {
		double spf = K_DESMOSOME;
		if (c1->agek > cont::THRESH_SP && c2->agek > cont::THRESH_SP) {
			spf = K_TOTAL;
		}
		return adhesion(c1, c2, spf);
	}
}

double CI_fix_mu_to_fix_mu::CIFuncCName(const Cell*const RESTRICT c1, const Cell*const RESTRICT c2) {
	//pair‚©‚Ç‚¤‚©‚ÍŠO‚Å”»’è‚µ‚æ‚¤
	if (is_near(c1, c2)) {
		return ljmain(c1, c2);
	}
	else {
		return adhesion(c1, c2, K_DESMOSOME);
	}
}

double CI_mu_to_memb::CIFuncCName(const Cell*const RESTRICT c1, const Cell*const RESTRICT c2) {
	//spring_force_to_memb‚ªÝ’è‚³‚ê‚Ä‚¢‚é•K—v‚ª‚ ‚é
	if (is_near(c1, c2)) {
		return ljmain(c1, c2);
	}
	else {
		return adhesion(c1, c2, c1->spring_force_to_memb);
	}
}

double CI_fix_to_memb::CIFuncCName(const Cell*const RESTRICT c1, const  Cell*const RESTRICT c2) {
	if (is_near(c1, c2)) {
		return ljmain(c1, c2);
	}
	else {
		const double distlj = sqrt(p_cell_dist_sq(c1, c2));
		const double LJ2 = distlj / (c1->radius + c2->radius);
		return -(Kspring / distlj)*(LJ2 - 1.0);
	}
}

double CI_memb_to_memb::CIFuncCName(const Cell*const RESTRICT c1, const  Cell*const RESTRICT c2) {
	
	

	using namespace cont;
	const double rad_sum = c1->radius+ c2->radius;
	const double rad_sum_sq = rad_sum*rad_sum;
	
	double cr_dist_sq;
	const double distSq = p_cell_dist_sq(c1, c2);
	if (distSq< (cr_dist_sq=rad_sum_sq*P_MEMB*P_MEMB)) {
		//assert(fabs(ljmain(me, oppo)) < 1000);
		const double LJ6 = POW3(cr_dist_sq) / POW3(distSq);
		//LJ6 = LJ6*LJ6*LJ6;
        double tmp = 4.0*eps_m*LJ6*(LJ6 - 1) / distSq;
        if (fabs(tmp)>100) {
            printf("menb too strong intr1,%d %d %lf %lf %lf\n",c1->get_index(),c2->get_index(), tmp, distSq, cr_dist_sq);
            assert(false);
            exit(1);
        }
		return tmp;
	}
	else if (distSq < rad_sum_sq) {
		const double distlj = sqrt(distSq);
		const double cr_dist = rad_sum*P_MEMB;
        double tmp = -(DER_DER_CONST / distlj) * (distlj / cr_dist - 1.0);
        if (fabs(tmp)>100) {
            printf("menb too strong intr2,%lf %lf %lf\n", tmp, distSq, cr_dist*cr_dist);
            assert(false);
            exit(1);
        }
		return tmp;
	}
	else {
		const double distlj = sqrt(distSq);
		const double lambda_dist = (1.0 + P_MEMB)*rad_sum;
		const double cr_dist = rad_sum*P_MEMB;
		const double LJ6 = POW6(cr_dist) / POW6(lambda_dist - distlj);
		//LJ6 = LJ6*LJ6;
		//LJ6 = LJ6*LJ6*LJ6;
        double tmp = -(DER_DER_CONST / rad_sum)*((1.0 - P_MEMB) / P_MEMB)
			- 4 * eps_m*(LJ6*(LJ6 - 1.0)) / ((lambda_dist - distlj)*distlj);
        if (fabs(tmp)>100) {
            printf("menb too strong intr3,%lf %lf %lf\n", tmp, distSq, cr_dist*cr_dist);
            assert(false);
            exit(1);
        }
        return tmp;
	}
}
double CI_memb_to_memb_diag::CIFuncCName(const Cell*const RESTRICT c1, const  Cell*const RESTRICT c2) {



    using namespace cont;
    const double rad_sum = c1->radius + c2->radius;
    const double rad_sum_sq = rad_sum*rad_sum;

    double cr_dist_sq;
    const double distSq = p_cell_dist_sq(c1, c2);
    if (distSq< (cr_dist_sq = rad_sum_sq*P_MEMB*P_MEMB)) {
        //assert(fabs(ljmain(me, oppo)) < 1000);
        const double LJ6 = POW3(cr_dist_sq) / POW3(distSq);
        //LJ6 = LJ6*LJ6*LJ6;
        double tmp = 4.0*eps_m*LJ6*(LJ6 - 1) / distSq;
       
        if (fabs(tmp)>100) {
            printf("menb too strong intr1,%d %d %lf %lf %lf\n", c1->get_index(), c2->get_index(), tmp, distSq, cr_dist_sq);
            assert(false);
            exit(1);
        }
        return tmp;
    }
    else if (distSq < rad_sum_sq) {
        const double distlj = sqrt(distSq);
        const double cr_dist = rad_sum*P_MEMB;
        double tmp = -(DER_DER_CONST / distlj) * (distlj / cr_dist - 1.0);
        if (fabs(tmp)>100) {
            printf("menb too strong intr2,%lf %lf %lf\n", tmp, distSq, cr_dist*cr_dist);
            assert(false);
            exit(1);
        }
        return tmp;
    }
    else {
        return 0;
    }
    
}

inline double CI_other::CIFuncCName(const Cell*const RESTRICT c1, const  Cell*const RESTRICT c2) {
	return ljmain(c1, c2);
}

inline double CI_pair::CIFuncCName(const Cell*const RESTRICT c1, const  Cell*const RESTRICT c2) {
	const double dist = sqrt(p_cell_dist_sq(c1, c2));
	const double force = Kspring_division*(dist - c1->spr_nat_len);
	return -force / dist;
}

inline void wall_interaction(Cell*const RESTRICT c) {
	const double distlj = 2.0*c->z();
	const double LJ6 = POW6(c->radius) / POW6(c->z());
	//LJ6 = LJ6*LJ6;
	//LJ6 = LJ6*LJ6*LJ6;
	const double ljm = 4.0*eps_m*LJ6*(LJ6 - 1.0) / (distlj*distlj);
	c->z += cont::DT_Cell* ljm*2.0*c->z();
}
////////////////////////////////////////////////////////////////////
inline void _memb_bend_calc1(Cell *const RESTRICT memb)
{
	const double dn = sqrt(p_cell_dist_sq(memb->md.memb_r, memb));

	memb->md.nv[0] = p_diff_x(memb->md.memb_r->x(), memb->x()) / dn;
	memb->md.nv[1] = p_diff_y(memb->md.memb_r->y(), memb->y()) / dn;
	memb->md.nv[2] = (memb->md.memb_r->z() - memb->z()) / dn;
	memb->md.dn = dn;

	const double dm = sqrt(p_cell_dist_sq(memb->md.memb_u, memb));
	memb->md.mv[0] = p_diff_x(memb->md.memb_u->x(), memb->x()) / dm;
	memb->md.mv[1] = p_diff_y(memb->md.memb_u->y(), memb->y()) / dm;
	memb->md.mv[2] = (memb->md.memb_u->z() - memb->z()) / dm;
	memb->md.dm = dm;

#ifdef DIAG_BEND
    const double dn_a = sqrt(p_cell_dist_sq(memb->md.memb_r->md.memb_b, memb)); //ur
    memb->md.nv_a[0] = p_diff_x(memb->md.memb_r->md.memb_b->x(), memb->x()) / dn_a;
    memb->md.nv_a[1] = p_diff_y(memb->md.memb_r->md.memb_b->y(), memb->y()) / dn_a;
    memb->md.nv_a[2] = (memb->md.memb_r->md.memb_b->z() - memb->z()) / dn_a;
    memb->md.dn_a = dn_a;

    const double dm_a = sqrt(p_cell_dist_sq(memb->md.memb_u->md.memb_r, memb)); //ur
    memb->md.mv_a[0] = p_diff_x(memb->md.memb_u->md.memb_r->x(), memb->x()) / dm_a;
    memb->md.mv_a[1] = p_diff_y(memb->md.memb_u->md.memb_r->y(), memb->y()) / dm_a;
    memb->md.mv_a[2] = (memb->md.memb_u->md.memb_r->z() - memb->z()) / dm_a;
    memb->md.dm_a = dm_a;
#endif
}

inline void _memb_bend_calc2(Cell *const RESTRICT memb)
{
	memb->md.ipn =
		memb->md.nv[0] * memb->md.memb_r->md.nv[0]
		+ memb->md.nv[1] * memb->md.memb_r->md.nv[1]
		+ memb->md.nv[2] * memb->md.memb_r->md.nv[2];

	memb->md.ipm =
		memb->md.mv[0] * memb->md.memb_u->md.mv[0]
		+ memb->md.mv[1] * memb->md.memb_u->md.mv[1]
		+ memb->md.mv[2] * memb->md.memb_u->md.mv[2];
#ifdef DIAG_BEND
    memb->md.ipn_a =
        memb->md.nv_a[0] * memb->md.memb_r->md.memb_b->md.nv_a[0]
        + memb->md.nv_a[1] * memb->md.memb_r->md.memb_b->md.nv_a[1]
        + memb->md.nv_a[2] * memb->md.memb_r->md.memb_b->md.nv_a[2];

    memb->md.ipm_a =
        memb->md.mv_a[0] * memb->md.memb_u->md.memb_r->md.mv_a[0]
        + memb->md.mv_a[1] * memb->md.memb_u->md.memb_r->md.mv_a[1]
        + memb->md.mv_a[2] * memb->md.memb_u->md.memb_r->md.mv_a[2];
#endif
}

inline void _memb_bend_calc3(Cell *const RESTRICT memb)
{




    auto& nvr = memb->md.memb_r->md.nv;
    auto& nvl = memb->md.memb_l->md.nv;
    auto& nvll = memb->md.memb_ll->md.nv;

    auto& ipnl = memb->md.memb_l->md.ipn;
    auto& ipnll = memb->md.memb_ll->md.ipn;


    auto& mvu = memb->md.memb_u->md.mv;
    auto& mvb = memb->md.memb_b->md.mv;
    auto& mvbb = memb->md.memb_bb->md.mv;

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

#define BCALC(suff,num)\
    (\
            -(1.0 - ipn##suff)*(nvr##suff[num] - ipn##suff*nv##suff[num]) / dn##suff\
            + (1.0 - ipnl##suff)*((nv##suff[num] - ipnl##suff*nvl##suff[num]) / dnl##suff - (nvl##suff[num] - ipnl##suff*nv##suff[num]) / dn##suff)\
            + (1.0 - ipnll##suff)*(nvll##suff[num] - ipnll##suff*nvl##suff[num]) / dnl##suff\
    \
            - (1.0 - ipm##suff)*(mvu##suff[num] - ipm##suff*mv##suff[num]) / dm##suff\
            + (1.0 - ipmb##suff)*((mv##suff[num] - ipmb##suff*mvb##suff[num]) / dmb##suff - (mvb##suff[num] - ipmb##suff*mv##suff[num]) / dm##suff)\
            + (1.0 - ipmbb##suff)*(mvbb##suff[num] - ipmbb##suff*mvb##suff[num]) / dmb##suff)

    memb->x +=cont::DT_Cell* KBEND*BCALC(,0);

    memb->y += cont::DT_Cell* KBEND*BCALC(,1);

    memb->z += cont::DT_Cell* KBEND*BCALC(,2);


#ifdef DIAG_BEND
    auto& nvr_a = memb->md.memb_r->md.memb_b->md.nv_a;
    auto& nvl_a = memb->md.memb_l->md.memb_u->md.nv_a;
    auto& nvll_a = memb->md.memb_l->md.memb_u->md.memb_l->md.memb_u->md.nv_a;
    auto& mvu_a = memb->md.memb_u->md.memb_r->md.mv_a;
    auto& mvb_a = memb->md.memb_b->md.memb_l->md.mv_a;
    auto& mvbb_a = memb->md.memb_b->md.memb_l->md.memb_b->md.memb_l->md.mv_a;



    auto& ipnl_a = memb->md.memb_l->md.memb_u->md.ipn_a;
    auto& ipnll_a = memb->md.memb_l->md.memb_u->md.memb_l->md.memb_u->md.ipn_a;

    auto& ipmb_a = memb->md.memb_b->md.memb_l->md.ipm_a;
    auto& ipmbb_a = memb->md.memb_b->md.memb_l->md.memb_b->md.memb_l->md.ipn_a;

    auto& dnl_a = memb->md.memb_l->md.memb_u->md.dn_a;
    auto& dmb_a = memb->md.memb_b->md.memb_l->md.dm_a;


    auto& dn_a = memb->md.dn_a;
    auto& dm_a = memb->md.dm_a;
    auto& nv_a = memb->md.nv_a;
    auto& mv_a = memb->md.mv_a;
    auto& ipn_a = memb->md.ipn_a;
    auto& ipm_a = memb->md.ipm_a;

    memb->x +=cont::DT_Cell* KBEND*BCALC(_a,0);

    memb->y += cont::DT_Cell* KBEND*BCALC(_a,1);

    memb->z += cont::DT_Cell* KBEND*BCALC(_a,2);
#endif



}

inline void memb_bend(CellManager& cman) {
	cman.memb_foreach_parallel_native([&](Cell*const RESTRICT c) {
		_memb_bend_calc1(c);
	});

	cman.memb_foreach_parallel_native([](Cell*const RESTRICT c) {
		_memb_bend_calc2(c);
	});

	cman.memb_foreach_parallel_native([](Cell*const RESTRICT c) {
		_memb_bend_calc3(c);
	});
}

//////////////////////////////////////////////////////////////////////////////////////

inline bool collide_with_wall(const Cell*const RESTRICT c) {
	return c->z() < c->radius;
}



inline void _MEMB_interaction(Cell*const RESTRICT memb) {
	if (collide_with_wall(memb)) {
		wall_interaction(memb);
	}

	size_t sz = memb->connected_cell.size();
	for (size_t i = 0; i < 4; i++) {
		if (memb->connected_cell[i]->state == MEMB&&no_double_count(memb, memb->connected_cell[i])) {
			cell_interaction_apply<CI_memb_to_memb>(memb, memb->connected_cell[i]);
		}
	}
    for (size_t i = 4; i < sz; i++) {
        if (memb->connected_cell[i]->state == MEMB&&no_double_count(memb, memb->connected_cell[i])) {
            cell_interaction_apply<CI_memb_to_memb_diag>(memb, memb->connected_cell[i]);
        }
    }
}

void _DER_interaction(Cell*const RESTRICT der) {
	if (collide_with_wall(der)) {
		wall_interaction(der);
	}
	der->connected_cell.foreach([&](Cell*const RESTRICT conn) {
		if (conn->state == DER && no_double_count(der, conn)) {
			cell_interaction_apply<CI_der_to_der>(der, conn);
		}
		else {
			cell_interaction_apply<CI_other>(der, conn);
		}
	});
}

void _AL_AIR_DE_interaction(Cell*const RESTRICT aad) {
	aad->connected_cell.foreach([&](Cell*const RESTRICT conn) {
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

inline void _FIX_interaction(Cell*const RESTRICT fix) {
    fix->connected_cell.foreach([&](Cell*const conn) { //cannot be restricted due to fix->pair
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
            if (fix->dermis()->memb_exist<cont::MEMB_ADHE_RANGE>(conn)) {
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


inline bool paired_with_fix(const Cell*const RESTRICT c) {
    if(c->pair==nullptr){
        return false;
    }
    return c->pair->state == FIX;
}
void _MUSUME_interaction(Cell*const RESTRICT& musume) {
	musume->spring_force_to_memb = 0;
	if (paired_with_fix(musume)) {
		musume->spring_force_to_memb = Kspring;
	}
	else if (musume->rest_div_times > 0) {
		musume->spring_force_to_memb = Kspring_d;
	}
    musume->connected_cell.foreach([&](Cell*const conn) { //cannot be restricted due to musume->pair
		if (conn->state == MEMB) {
            if (musume->dermis()->memb_exist<cont::MEMB_ADHE_RANGE>(conn)) {
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


inline void _pair_interaction(Cell*const RESTRICT paired) {
	if (no_double_count(paired, paired->pair)) {
		cell_interaction_apply<CI_pair>(paired, paired->pair);
	}
}


//////////////////////////////////////////////////////////////////////////////////////



void cell_interaction(CellManager & cman)
{
    memb_bend(cman);
   // printf("check...\n");
		cman.memb_foreach_parallel_native([](Cell*const RESTRICT c) {
			_MEMB_interaction(c);
		});

		cman.non_memb_foreach_parallel_native([](Cell*const RESTRICT c) {
			switch (c->state) {
				/*
				case MEMB:
				_MEMB_interaction(c);
				break;
				*/
			case DER:
				_DER_interaction(c);
				break;
			case FIX:
				_FIX_interaction(c);
				break;
			case MUSUME:
				_MUSUME_interaction(c);
				break;
			case ALIVE:case AIR:case DEAD:
				_AL_AIR_DE_interaction(c);
				break;
			default:
				break;
			}
			if (c->pair != nullptr)_pair_interaction(c);
		});
		
		
		

	
}
