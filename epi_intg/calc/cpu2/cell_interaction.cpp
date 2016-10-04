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
#include <cmath>

#define POW6(x) ((x)*(x)*(x)*(x)*(x)*(x))
#define POW3(x) ((x)*(x)*(x))
#define CIFuncCName Coef
#define CIFuncDecl(name) struct name{static real CIFuncCName(const Cell*const RESTRICT c1,const Cell*const RESTRICT c2);}

template<class Fn>
void cell_interaction_apply(Cell*const RESTRICT c1, Cell*const RESTRICT c2) {
	const real coef = Fn::CIFuncCName(c1, c2);
    if (fabs(coef) > 1000) {
    	throw std::logic_error("Too strong interaction between following two cells:\n"_s
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
	return 4.0*pm->eps_m*LJ6*(LJ6 - 1) / (dist*dist2) + pm->para_ljp2;
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
	return 4.0*pm->eps_m*LJ6*(LJ6 - 1) / (dist*dist) + pm->para_ljp2;
}

inline real ljmain(const Cell*const RESTRICT c1, const  Cell*const RESTRICT c2) {

	const real distSq = p_cell_dist_sq(c1, c2);
	const real rad_sum = c1->radius + c2->radius;
	const real LJ6 = POW6(rad_sum) / POW3(distSq);
	//LJ6 = LJ6*LJ6*LJ6;
	return 4.0*pm->eps_m*LJ6*(LJ6 - 1) / distSq;

}

inline real adhesion(const Cell*const RESTRICT c1, const  Cell*const RESTRICT c2,const real spring_const) {
	const real distlj = sqrt(p_cell_dist_sq(c1, c2));
	const real rad_sum = c1->radius + c2->radius;
	//real LJ2_m1;

	const real LJ2_m1 = (distlj / rad_sum) - 1;
	return (LJ2_m1 + 1 > pm->LJ_THRESH ?
		0.0 :
		-(spring_const / distlj)
		* LJ2_m1*(1 - LJ2_m1*LJ2_m1 / ((pm->LJ_THRESH - 1.0)*(pm->LJ_THRESH - 1.0))));
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
		return -(pm->Kspring / distlj)*(LJ2 - 1.0);
	}
}

real CI_memb_to_memb::CIFuncCName(const Cell*const RESTRICT c1, const  Cell*const RESTRICT c2) {
	const real rad_sum = c1->radius+ c2->radius;
	const real rad_sum_sq = rad_sum*rad_sum;

	real cr_dist_sq;
	const real distSq = p_cell_dist_sq(c1, c2);
	if (distSq< (cr_dist_sq=rad_sum_sq*pm->P_MEMB*pm->P_MEMB)) {
		const real LJ6 = POW3(cr_dist_sq) / POW3(distSq);
        return 4.0*pm->eps_m*LJ6*(LJ6 - 1) / distSq;
	}
	else if (distSq < rad_sum_sq) {
		const real distlj = sqrt(distSq);
		const real cr_dist = rad_sum*pm->P_MEMB;
        return -(pm->DER_DER_CONST / distlj) * (distlj / cr_dist - 1.0);
	}
	else {
		const real distlj = sqrt(distSq);
		const real lambda_dist = (1.0 + pm->P_MEMB)*rad_sum;
		const real cr_dist = rad_sum*pm->P_MEMB;
		const real LJ6 = POW6(cr_dist) / POW6(lambda_dist - distlj);
        return -(pm->DER_DER_CONST / rad_sum)*((1.0 - pm->P_MEMB) / pm->P_MEMB)
			- 4 * pm->eps_m*(LJ6*(LJ6 - 1.0)) / ((lambda_dist - distlj)*distlj);
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
	const real distlj = 2.0*c->z();
	const real LJ6 = POW6(c->radius) / POW6(c->z());
	//LJ6 = LJ6*LJ6;
	//LJ6 = LJ6*LJ6*LJ6;
	const real ljm = 4.0*pm->eps_m*LJ6*(LJ6 - 1.0) / (distlj*distlj);
	c->z += pm->DT_Cell* ljm*2.0*c->z();
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

