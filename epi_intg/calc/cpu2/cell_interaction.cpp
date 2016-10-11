/*
 * cell_interaction.cpp
 *
 *  Created on: 2016/10/04
 *      Author: yasu7890v
 */

#include "cell_interaction.h"
#include "../../utils.h"
#include "../../define.h"
#include "../../global.h"
#include "cell.h"
#include "CellManager.h"
#include <cmath>

#define POW6(x) ((x)*(x)*(x)*(x)*(x)*(x))
#define POW3(x) ((x)*(x)*(x))
#define CIFuncCName Coef
#define CIFuncDecl(name) struct name{static real CIFuncCName(const Cell*const RESTRICT c1,const Cell*const RESTRICT c2);}

template<class Fn>
void cell_interaction_apply(Cell*const RESTRICT c1, Cell*const RESTRICT c2) {
	const real coef = Fn::CIFuncCName(c1, c2);
    if (fabs(coef) > 1000) {
    	throw std::logic_error("Too strong interaction(coef="_s+std::to_string(coef)+") between following two cells:\n"_s
    			+"1.\n"+c1->cell_info_str()+"\n2.\n"+c2->cell_info_str());
    }
	const real dumx = pm->DT_Cell*coef*p_diff_x(c1->x(), c2->x());
	const real dumy = pm->DT_Cell*coef*p_diff_y(c1->y(), c2->y());
	const real dumz = pm->DT_Cell*coef*(c1->z() - c2->z());
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
	return sqrt(p_cell_dist_sq(der1, der2))+ pm->delta_R < der1->radius + der2->radius;
}

inline real ljmain_der_near(const Cell*const RESTRICT der1, const Cell*const RESTRICT der2) {
	const real dist = sqrt(p_cell_dist_sq(der1, der2));
	const real dist2 = dist + pm->delta_R;
	const real rad_sum = der1->radius + der2->radius;
	const real LJ6 = POW6(rad_sum) / POW6(dist2);
	//LJ6 = LJ6*LJ6*LJ6;
	return real(4.0)*pm->eps_m*LJ6*(LJ6 - real(1.0)) / (dist*dist2) + pm->para_ljp2;
}

inline bool is_near(const Cell*const RESTRICT c1, const Cell*const RESTRICT c2) {
	//real rad_sum = c1->radius + c2->radius;
	return p_cell_dist_sq(c1, c2) < (c1->radius + c2->radius)*(c1->radius + c2->radius);
}

inline real ljmain_der_far(const Cell*const RESTRICT der1, const Cell*const RESTRICT der2) {
	const real dist = sqrt(p_cell_dist_sq(der1, der2));
	const real dist2 = dist + pm->delta_R;
	const real rad_sum = der1->radius + der2->radius;
	const real LJ6 = POW6(rad_sum) / POW6(dist2);
	return real(4.0)*pm->eps_m*LJ6*(LJ6 - real(1.0)) / (dist*dist) + pm->para_ljp2;
}

inline real ljmain(const Cell*const RESTRICT c1, const  Cell*const RESTRICT c2) {

	const real distSq = p_cell_dist_sq(c1, c2);
	const real rad_sum = c1->radius + c2->radius;
	const real LJ6 = POW6(rad_sum) / POW3(distSq);
	//LJ6 = LJ6*LJ6*LJ6;
	return real(4.0)*pm->eps_m*LJ6*(LJ6 - real(1.0)) / distSq;

}

inline real adhesion(const Cell*const RESTRICT c1, const  Cell*const RESTRICT c2,const real spring_const) {
	const real distlj = sqrt(p_cell_dist_sq(c1, c2));
	const real rad_sum = c1->radius + c2->radius;
	//real LJ2_m1;

	const real LJ2_m1 = (distlj / rad_sum) - real(1.0);
	return (LJ2_m1 + real(1.0) > pm->LJ_THRESH ?
        real(0.0) :
		-(spring_const / distlj)
		* LJ2_m1*(real(1.0) - LJ2_m1*LJ2_m1 / ((pm->LJ_THRESH - real(1.0))*(pm->LJ_THRESH - real(1.0)))));
}

/////////////////////////////////////////////////////////////////////////


real CI_der_to_der::CIFuncCName(const Cell*const RESTRICT c1, const Cell*const RESTRICT c2) {
	if (is_der_near(c1, c2)) {
		return ljmain_der_near(c1, c2);
	}
	else if (is_near(c1, c2)) {
		return pm->para_ljp2;
	}
	else {
		return ljmain_der_far(c1, c2);
	}
}

real CI_al_air_de_to_al_air_de_fix_mu::CIFuncCName(const Cell*const RESTRICT c1, const Cell*const RESTRICT c2) {
	if (is_near(c1, c2)) {
		return ljmain(c1, c2);
	}
	else {
		real spf = pm->K_DESMOSOME;
		if (c1->agek > pm->THRESH_SP && c2->agek > pm->THRESH_SP) {
			spf = pm->K_TOTAL;
		}
		return adhesion(c1, c2, spf);
	}
}

real CI_fix_mu_to_fix_mu::CIFuncCName(const Cell*const RESTRICT c1, const Cell*const RESTRICT c2) {
	//pairかどうかは外で判定しよう
	if (is_near(c1, c2)) {
		return ljmain(c1, c2);
	}
	else {
		return adhesion(c1, c2, pm->K_DESMOSOME);
	}
}

real CI_mu_to_memb::CIFuncCName(const Cell*const RESTRICT c1, const Cell*const RESTRICT c2) {
	//spring_force_to_membが設定されている必要がある
	if (is_near(c1, c2)) {
		return ljmain(c1, c2);
	}
	else {
		return adhesion(c1, c2, c1->spring_force_to_memb);
	}
}

real CI_fix_to_memb::CIFuncCName(const Cell*const RESTRICT c1, const  Cell*const RESTRICT c2) {
	if (is_near(c1, c2)) {
		return ljmain(c1, c2);
	}
	else {
		const real distlj = sqrt(p_cell_dist_sq(c1, c2));
		const real LJ2 = distlj / (c1->radius + c2->radius);
		return -(pm->Kspring / distlj)*(LJ2 - real(1.0));
	}
}

real CI_memb_to_memb::CIFuncCName(const Cell*const RESTRICT c1, const  Cell*const RESTRICT c2) {
	const real rad_sum = c1->radius+ c2->radius;
	const real rad_sum_sq = rad_sum*rad_sum;

	real cr_dist_sq;
	const real distSq = p_cell_dist_sq(c1, c2);
	if (distSq< (cr_dist_sq=rad_sum_sq*pm->P_MEMB*pm->P_MEMB)) {
		const real LJ6 = POW3(cr_dist_sq) / POW3(distSq);
        return real(4.0)*pm->eps_m*LJ6*(LJ6 - real(1.0)) / distSq;
	}
	else if (distSq < rad_sum_sq) {
		const real distlj = sqrt(distSq);
		const real cr_dist = rad_sum*pm->P_MEMB;
        return -(pm->DER_DER_CONST / distlj) * (distlj / cr_dist - real(1.0));
	}
	else {
		const real distlj = sqrt(distSq);
		const real lambda_dist = (real(1.0) + pm->P_MEMB)*rad_sum;
		const real cr_dist = rad_sum*pm->P_MEMB;
		const real LJ6 = POW6(cr_dist) / POW6(lambda_dist - distlj);
        return -(pm->DER_DER_CONST / rad_sum)*((real(1.0) - pm->P_MEMB) / pm->P_MEMB)
			- real(4.0) * pm->eps_m*(LJ6*(LJ6 - real(1.0))) / ((lambda_dist - distlj)*distlj);
	}
}


inline real CI_other::CIFuncCName(const Cell*const RESTRICT c1, const  Cell*const RESTRICT c2) {
	return ljmain(c1, c2);
}

inline real CI_pair::CIFuncCName(const Cell*const RESTRICT c1, const  Cell*const RESTRICT c2) {
	const real dist = sqrt(p_cell_dist_sq(c1, c2));
	const real force = pm->Kspring_division*(dist - c1->spr_nat_len);
	return -force / dist;
}

inline void wall_interaction(Cell*const RESTRICT c) {
	const real distlj = real(2.0)*c->z();
	const real LJ6 = POW6(c->radius) / POW6(c->z());
	//LJ6 = LJ6*LJ6;
	//LJ6 = LJ6*LJ6*LJ6;
	const real ljm = real(4.0)*pm->eps_m*LJ6*(LJ6 - real(1.0)) / (distlj*distlj);
	c->z += pm->DT_Cell* ljm*real(2.0)*c->z();
}
////////////////////////////////////////////////////////////////////

inline void _memb_bend_calc1(Cell *const RESTRICT memb)
{
	const real dn = sqrt(p_cell_dist_sq(memb->md.memb[0], memb));

	memb->md.nv[0] = p_diff_x(memb->md.memb[0]->x(), memb->x()) / dn;
	memb->md.nv[1] = p_diff_y(memb->md.memb[0]->y(), memb->y()) / dn;
	memb->md.nv[2] = (memb->md.memb[0]->z() - memb->z()) / dn;
	memb->md.dn = dn;

	const real dm = sqrt(p_cell_dist_sq(memb->md.memb[1], memb));
	memb->md.mv[0] = p_diff_x(memb->md.memb[1]->x(), memb->x()) / dm;
	memb->md.mv[1] = p_diff_y(memb->md.memb[1]->y(), memb->y()) / dm;
	memb->md.mv[2] = (memb->md.memb[1]->z() - memb->z()) / dm;
	memb->md.dm = dm;

}

inline void _memb_bend_calc2(Cell *const RESTRICT memb)
{
	memb->md.ipn =
		memb->md.nv[0] * memb->md.memb[0]->md.nv[0]
		+ memb->md.nv[1] * memb->md.memb[0]->md.nv[1]
		+ memb->md.nv[2] * memb->md.memb[0]->md.nv[2];

	memb->md.ipm =
		memb->md.mv[0] * memb->md.memb[1]->md.mv[0]
		+ memb->md.mv[1] * memb->md.memb[1]->md.mv[1]
		+ memb->md.mv[2] * memb->md.memb[1]->md.mv[2];
}

struct MEMB_bend_diff_factor_orig {
    real operator()(real cdot) {
        return (real(1.0) - cdot);
    }
};

struct MEMB_bend_diff_factor_new {
    real operator()(real cdot) {
        if (cdot >= real(0.0)) {
            return (real(1.0) - cdot);
        }
        else {
            return real(1.0) / ((real(1.0) + cdot)*(real(1.0) + cdot));
        }
    }
};
template<typename Fn>
inline void _memb_bend_calc3(Cell *const RESTRICT memb,Fn _memb_bend_diff_factor)
{




    auto& nvr = memb->md.memb[0]->md.nv;
    auto& nvl = memb->md.memb[2]->md.nv;
    auto& nvll = memb->md.memb[2]->md.memb[2]->md.nv;

    auto& ipnl = memb->md.memb[2]->md.ipn;
    auto& ipnll = memb->md.memb[2]->md.memb[2]->md.ipn;


    auto& mvu = memb->md.memb[1]->md.mv;
    auto& mvb = memb->md.memb[3]->md.mv;
    auto& mvbb = memb->md.memb[3]->md.memb[3]->md.mv;

    auto& ipmb = memb->md.memb[3]->md.ipm;
    auto& ipmbb = memb->md.memb[3]->md.memb[3]->md.ipn;



    auto& dnl = memb->md.memb[0]->md.dn;
    auto& dmb = memb->md.memb[3]->md.dm;



    auto& dn = memb->md.dn;
    auto& dm = memb->md.dm;
    auto& nv = memb->md.nv;
    auto& mv = memb->md.mv;
    auto& ipn = memb->md.ipn;
    auto& ipm = memb->md.ipm;

#define BCALC(suff,num)\
    (\
            -(_memb_bend_diff_factor(ipn##suff))*(nvr##suff[num] - ipn##suff*nv##suff[num]) / dn##suff\
            + _memb_bend_diff_factor(ipnl##suff)*((nv##suff[num] - ipnl##suff*nvl##suff[num]) / dnl##suff - (nvl##suff[num] - ipnl##suff*nv##suff[num]) / dn##suff)\
            + _memb_bend_diff_factor(ipnll##suff)*(nvll##suff[num] - ipnll##suff*nvl##suff[num]) / dnl##suff\
    \
            - _memb_bend_diff_factor(ipm##suff)*(mvu##suff[num] - ipm##suff*mv##suff[num]) / dm##suff\
            + _memb_bend_diff_factor(ipmb##suff)*((mv##suff[num] - ipmb##suff*mvb##suff[num]) / dmb##suff - (mvb##suff[num] - ipmb##suff*mv##suff[num]) / dm##suff)\
            + _memb_bend_diff_factor(ipmbb##suff)*(mvbb##suff[num] - ipmbb##suff*mvb##suff[num]) / dmb##suff)

    memb->x += pm->DT_Cell* pm->KBEND*BCALC(, 0);

    memb->y += pm->DT_Cell* pm->KBEND*BCALC(, 1);

    memb->z += pm->DT_Cell* pm->KBEND*BCALC(, 2);



}

Vec<3,real> tri_normal(const Cell* const RESTRICT c1, const Cell* const RESTRICT c2, const Cell* const RESTRICT c3) {
    Vec<3, real> c1v = c1->pos_as_vector();
    Vec<3, real> c2v = c2->pos_as_vector();
    Vec<3, real> c3v = c3->pos_as_vector();
    return Vec<3, real>::cross(c1v - c3v,c2v - c3v).normalize();
}

inline Vec<3, real> tri_rho(const Vec<3, real>& n1, const Vec<3, real>& b1) {
    return b1 - (Vec<3, real>::dot(n1,b1)*n1);
}

template<class Fn>
Vec<3, real>  tri_bend_1(const Cell* const RESTRICT r1, const Cell* const RESTRICT a1, const Cell* const RESTRICT b1, const Cell* const RESTRICT c1,Fn _memb_bend_diff_factor) {
    const Vec<3, real>  r1v = r1->pos_as_vector(), a1v = a1->pos_as_vector(), b1v = b1->pos_as_vector(), c1v = c1->pos_as_vector();
    const Vec<3, real> sub1 = vpsub(a1v, b1v);
    Vec<3, real> norm1 = Vec<3, real>::cross(vpsub(r1v, b1v), sub1);

    real norm1_norm = norm1.normalize_with_norm();
    Vec<3, real> norm2 = Vec<3, real>::cross(sub1, vpsub(c1v, b1v)).normalize();
    return _memb_bend_diff_factor(Vec<3, real>::dot(norm1, norm2))*(Vec<3, real>::cross(sub1, tri_rho(norm1, norm2)) / norm1_norm);
}

template<class Fn>
Vec<3, real>  tri_bend_2(const Cell* const RESTRICT r1, const Cell* const RESTRICT a1, const Cell* const RESTRICT b1, const Cell* const RESTRICT c1, Fn _memb_bend_diff_factor) {
    const Vec<3, real>  r1v = r1->pos_as_vector(), a1v = a1->pos_as_vector(), b1v = b1->pos_as_vector(), c1v = c1->pos_as_vector();
    const Vec<3, real> sub1 = vpsub(b1v, c1v), sub2 = vpsub(a1v, b1v);
    Vec<3, real> norm1 = Vec<3, real>::cross(vpsub(r1v, c1v), sub1);
    Vec<3, real> norm2 = Vec<3, real>::cross(vpsub(r1v, b1v), sub2);

    real norm1_norm = norm1.normalize_with_norm(), norm2_norm = norm2.normalize_with_norm();
    return _memb_bend_diff_factor(Vec<3, real>::dot(norm1, norm2))*(Vec<3, real>::cross(sub1, tri_rho(norm1, norm2)) / norm1_norm + Vec<3, real>::cross(sub2, tri_rho(norm2, norm1)) / norm2_norm);
}

template<class Fn>
void tri_memb_bend(Cell* const RESTRICT memb, Fn _memb_bend_diff_factor) {
    Vec<3, real> dr = { 0.0,0.0,0.0 };
    const size_t msz = memb->md.memb.size();
    for (int i = 0; i<msz; i++) {
        dr += tri_bend_1(memb, memb->md.memb[i], memb->md.memb[(i + 1) % msz], memb->md.memb[(i + 1) % msz]->md.memb[i], _memb_bend_diff_factor);

    }
    for (int i = 0; i<msz; i++) {
        dr += tri_bend_2(memb, memb->md.memb[i], memb->md.memb[(i + 1) % msz], memb->md.memb[(i + 2) % msz], _memb_bend_diff_factor);
    }

    memb->x += pm->DT_Cell*pm->KBEND*dr[0];
    memb->y += pm->DT_Cell*pm->KBEND*dr[1];
    memb->z += pm->DT_Cell*pm->KBEND*dr[2];

}

template<class Fn>
inline void memb_bend(CellManager& cman, Fn _memb_bend_diff_factor) {
    if (pm->USE_TRI_MEMB) {
        cman.memb_foreach_parallel_native([&](Cell*const RESTRICT c) {
            tri_memb_bend(c,_memb_bend_diff_factor);
        });
    }
    else {
        cman.memb_foreach_parallel_native([&](Cell*const RESTRICT c) {
            _memb_bend_calc1(c);
        });

        cman.memb_foreach_parallel_native([](Cell*const RESTRICT c) {
            _memb_bend_calc2(c);
        });

        cman.memb_foreach_parallel_native([&](Cell*const RESTRICT c) {
            _memb_bend_calc3(c, _memb_bend_diff_factor);
        });
    }
}

inline bool collide_with_wall(const Cell*const RESTRICT c) {
    return c->z() < c->radius;
}



inline void _MEMB_interaction(Cell*const RESTRICT memb) {
    if (collide_with_wall(memb)) {
        wall_interaction(memb);
    }

    for (size_t i = 0; i < pm->MEMB_ADJ_CONN_NUM; i++) {
        if (memb->connected_cell[i]->state == MEMB&&no_double_count(memb, memb->connected_cell[i])) {
            cell_interaction_apply<CI_memb_to_memb>(memb, memb->connected_cell[i]);
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
            if (fix->pair != conn&&no_double_count(fix, conn)) {
                cell_interaction_apply<CI_fix_mu_to_fix_mu>(fix, conn);
            }
            break;
        case MUSUME:
            if (fix->pair != conn)cell_interaction_apply<CI_fix_mu_to_fix_mu>(fix, conn);
            break;
        case MEMB:
            if (fix->dermis()==conn) {
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
    if (c->pair == nullptr) {
        return false;
    }
    return c->pair->state == FIX;
}
void _MUSUME_interaction(Cell*const RESTRICT& musume) {
    musume->spring_force_to_memb = 0;
    if (paired_with_fix(musume)) {
        musume->spring_force_to_memb = pm->Kspring;
    }
    else if (musume->rest_div_times > 0) {
        musume->spring_force_to_memb = pm->Kspring_d;
    }
    musume->connected_cell.foreach([&](Cell*const conn) { //cannot be restricted due to musume->pair
        if (conn->state == MEMB) {
            if (musume->dermis()==conn) {
                cell_interaction_apply<CI_mu_to_memb>(musume, conn);
            }
            else {
                cell_interaction_apply<CI_other>(musume, conn);
            }
        }
        else if (conn->state == MUSUME &&musume->pair != conn&& no_double_count(musume, conn)) {
            cell_interaction_apply<CI_fix_mu_to_fix_mu>(musume, conn);
        }
    });
}


inline void _pair_interaction(Cell*const RESTRICT paired) {
    if (no_double_count(paired, paired->pair)) {
        cell_interaction_apply<CI_pair>(paired, paired->pair);
    }
}

void cell_interaction(CellManager & cman)
{

    if (pm->NEW_BEND_POT) {
        memb_bend(cman,MEMB_bend_diff_factor_orig());
    }
    else {
        memb_bend(cman, MEMB_bend_diff_factor_new());
    }

    cman.memb_foreach_parallel_native([](Cell*const RESTRICT c) {
        _MEMB_interaction(c);
    });

    cman.non_memb_foreach_parallel_native([](Cell*const RESTRICT c) {
        switch (c->state) {
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
