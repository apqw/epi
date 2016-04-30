#include "cell.h"
#include "funcs.h"
#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>
#include "static_if.hpp"


void calc_normal(const CellData& me_data,const CellData& dermis_data, double *nx, double *ny, double *nz) {
	//assert(me_data.dermis != nullptr);
	double nx1 = perio_diff_x(me_data.x, dermis_data.x);
	double ny1 = perio_diff_y(me_data.y, dermis_data.y);
	double nz1 = me_data.z - dermis_data.z;
	double norm = sqrt(normSq(nx1, ny1, nz1));

	*nx = nx1 / norm;
	*ny = ny1 / norm;
	*nz = nz1 / norm;

}

void Cell::div_direction(const CellData& me_data, const CellData& dermis_data, double *nx, double *ny, double *nz) {
	double rand_arg, rand_theta, rand_phi, cr1, sr1, cr2, sr2;
	double randnx, randny, randnz, sum;
	double nx1, ny1, nz1;
	double lx, ly, lz;

	calc_normal(me_data,dermis_data,&lx, &ly, &lz);
	do {
		rand_theta = M_PI*genrand_real();
		cr1 = cos(rand_theta);
		sr1 = sin(rand_theta);
		rand_phi = 2 * M_PI*genrand_real();
		cr2 = cos(rand_phi);
		sr2 = sin(rand_phi);

		randnx = sr1 * cr2;
		randny = sr1 * sr2;
		randnz = cr1;

		// take the outer product of the normal l and a random unit vector
		nx1 = ly * randnz - lz * randny;
		ny1 = lz * randnx - lx * randnz;
		nz1 = lx * randny - ly * randnx;
	} while ((sum = nx1*nx1 + ny1*ny1 + nz1*nz1) < 1.0e-14);
	sum = sqrt(sum);
	*nx = nx1 / sum;
	*ny = ny1 / sum;
	*nz = nz1 / sum;
}
void interac_wall(CellPtr& me) {
	double distlj, LJ1, dum_lj6, dum_ljmain,zz;
	zz = me->old_data.z;
	distlj = 2 * zz;
	//double d_inv = 1 / distlj;
	LJ1 = 2 * me->old_data.radius/distlj;
	double lj1sq = LJ1*LJ1;
	dum_lj6 = lj1sq*lj1sq*lj1sq;
	dum_ljmain = 4.0*eps_m * dum_lj6 * (dum_lj6 - 1.0) /(distlj*distlj);

	me->data.z = zz + DT_Cell*(4.0*eps_m * dum_lj6 * (dum_lj6 - 1.0) / (distlj*distlj)) * 2 * zz;
}

/*
	This ptr is transient. (this points something in std::vector)
	Do not keep long.
*/
CellPtr find_dermis(const CellPtr& me) {
	CellPtr* cptr_ptr = nullptr;
	double d1Sq = LX*LX;
	double distSq;
	auto& conn_memb_cell_set = me->old_data.connected_cell[MEMB];
	for (auto& ccell : conn_memb_cell_set) {
		distSq = dist3Sq(me->old_data.x, me->old_data.y, me->old_data.z,
			ccell->old_data.x, ccell->old_data.y, ccell->old_data.z);
		if (distSq < d1Sq) {
			d1Sq = distSq;
			cptr_ptr = &ccell;
			//if you use CellPtr directly,copy occurs
		}
	}
	if (cptr_ptr != nullptr) {
		return *cptr_ptr;
	}
	else {
		return nullptr;
	}
	//do not treat ptr directly, 
}


double k_lipid_release(double u, double age)
{
	return 0.25*lipid_rel*(1 + tanh((u - ubar) / 0.01))*(1 + tanh((age - THRESH_SP) / delta_lipid));
}

double k_lipid(double u, double age)
{
	return 0.25*lipid*(1 + tanh((ubar - u) / 0.01))*(1 + tanh((age - THRESH_SP) / delta_lipid));
}

double memb_memb_coef::operator()(const CellPtr& me, const CellPtr& oppo, double diffx, double diffy, double diffz) {
	/* core of heavy part */
	//		std::cout << "hello" << std::endl;
	double distSq = normSq(diffx, diffy, diffz);
	double rad_sum = me->old_data.radius + oppo->old_data.radius;
	double rad_sumSq = rad_sum*rad_sum;
	double cr_dist = rad_sum*P_MEMB; double cr_distSq = rad_sumSq*P_MEMB*P_MEMB;//constexpr

	if (distSq < cr_distSq) {
		double LJ6 = cr_distSq*cr_distSq*cr_distSq / (distSq*distSq*distSq);
		return 4 * eps_m * LJ6*(LJ6 - 1.0) / distSq;
	}
	else {
		double distlj = sqrt(distSq);
		if (distSq < rad_sumSq) {

			return -(DER_DER_CONST / distlj) * (distlj / cr_dist - 1.0);
		}
		else {
			double lambda_dist = (1.0 + P_MEMB)*rad_sum;
			double dist_delta = lambda_dist - distlj; double dist_deltaSq = dist_delta*dist_delta;
			double dum_lj6 = cr_distSq*cr_distSq*cr_distSq / (dist_deltaSq*dist_deltaSq*dist_deltaSq);

			return -(DER_DER_CONST / rad_sum)*((1.0 - P_MEMB) / P_MEMB)
				- 4 * eps_m*(dum_lj6*(dum_lj6 - 1.0)) / (dist_delta*distlj);
		}
	}
}

double null_coef::operator()(const CellPtr& me, const CellPtr& oppo, double diffx, double diffy, double diffz) {
	return 0;
}

double common_common_coef::operator()(const CellPtr& me, const CellPtr& oppo, double diffx, double diffy, double diffz) {
	double distSq = normSq(diffx, diffy, diffz);
	double LJ1Sq = (me->old_data.radius + oppo->old_data.radius) / distSq;
	double LJ6 = LJ1Sq*LJ1Sq*LJ1Sq;
	return 4.0*eps_m * LJ6 * (LJ6 - 1.0) / distSq;
}

double der_der_coef::operator()(const CellPtr& me, const CellPtr& oppo, double diffx, double diffy, double diffz) {
	double rad_sum = me->old_data.radius + oppo->old_data.radius;
	double rad_sum_d = rad_sum - delta_R;

	double rad_sumSq = rad_sum*rad_sum;
	double rad_sum_d_Sq = rad_sum_d*rad_sum_d;
	double distSq = normSq(diffx, diffy, diffz);
	if (distSq < rad_sum_d_Sq) {
		double dist = sqrt(distSq);
		double dist_d = dist + delta_R;
		double dist_d_Sq = dist_d*dist_d;
		double LJ1Sq = rad_sumSq / dist_d_Sq;
		double LJ6 = LJ1Sq*LJ1Sq*LJ1Sq;
		return 4.0*eps_m * LJ6 * (LJ6 - 1.0) / (dist_d*dist) + para_ljp2;
	}
	else if (distSq<rad_sumSq) {
		return para_ljp2;
	}
	else {
		double LJ1Sq = rad_sumSq / distSq;
		double LJ6 = LJ1Sq*LJ1Sq*LJ1Sq;
		return 4.0*eps_m * LJ6 * (LJ6 - 1.0) / distSq + para_ljp2;
	}
}

double supra_others_coef::operator()(const CellPtr& me, const CellPtr& oppo, double diffx, double diffy, double diffz) {
	double rad_sum = me->old_data.radius + oppo->old_data.radius;
	double rad_sumSq = rad_sum*rad_sum;
	double distSq = normSq(diffx, diffy, diffz);
	if (distSq < rad_sumSq) {
		double LJ1Sq = rad_sumSq / distSq;
		double LJ6 = LJ1Sq*LJ1Sq*LJ1Sq;
		return 4.0*eps_m * LJ6 * (LJ6 - 1.0) / distSq;
	}
	else {
		return adhesion(sqrt(distSq), rad_sum,
			me->old_data.agek>THRESH_SP && oppo->old_data.agek > THRESH_SP ?
			K_TOTAL :
			K_DESMOSOME);
	}
}

double stem_stem_coef::operator()(const CellPtr& me, const CellPtr& oppo, double diffx, double diffy, double diffz) {
	if (me->old_data.div_pair_cell == oppo) {
		return 0;
	}
	else {
		double rad_sum = me->old_data.radius + oppo->old_data.radius;

		double rad_sumSq = rad_sum*rad_sum;

		double distSq = normSq(diffx, diffy, diffz);
		if (distSq < rad_sumSq) {
			double LJ1Sq = rad_sumSq / distSq;
			double LJ6 = LJ1Sq*LJ1Sq*LJ1Sq;
			return 4.0*eps_m * LJ6 * (LJ6 - 1.0) / distSq;
		}
		else {
			return adhesion(sqrt(distSq), rad_sum, K_DESMOSOME);
		}
	}
}

double pair_coef::operator()(const CellPtr& me, const CellPtr& oppo, double diffx, double diffy, double diffz) {
	double dist = sqrt(normSq(diffx, diffy, diffz));
	double force = Kspring_division*(dist - me->old_data.pair_spr_nat_len);
	return force / dist;
}