#include "utils.h"
#include "CellManager.h"
#include "cell_interaction.h"
#include <cstdio>
#include <device_functions.h>
#define POW6(x) ((x)*(x)*(x)*(x)*(x)*(x))
#define POW3(x) ((x)*(x)*(x))

/** LJƒ|ƒeƒ“ƒVƒƒƒ‹‚ÌŒW” */
#define eps_m (0.01f)

/*
MEMB_calc
*/
#define P_MEMB (1.0f/COMPRESS_FACTOR)

/** L‚Ñ’e«ŒW” */
#define DER_DER_CONST (0.05f)
#define K_TOTAL (3.0f)
#define K_DESMOSOME_RATIO (0.01f)
#define K_DESMOSOME (3.0f*0.01f)
#define Kspring (25.0f)
#define Kspring_d (5.0f)

/** ‹È‚°’e«ŒW” */
#define KBEND (0.25f) //0.5->0.25

/*
DER_calc
*/
#define delta_R (0.4f*R_der)

#define para_ljp2 (0.005f)
#define Kspring_division (5.0f)

__device__ inline real4 normal_with_len_inv(const real4 a, const real4 b){
	const real diffx = p_diff_x(a.x, b.x);
	const real diffy = p_diff_y(a.y, b.y);
	const real diffz = a.z - b.z;
	const real len_inv = rsqrtr(diffx*diffx + diffy*diffy + diffz*diffz);
	return make_real4(diffx*len_inv, diffy*len_inv, diffz*len_inv, len_inv);
}

template<unsigned int AXIS>
__device__ inline real bend_calc_1(const real dot_rme, const  real dot_l, const real dot_ll, const real4 n_r, const real4 n_rme, const  real4 n_l, const  real4 n_ll){
	return 0.0f;
}//dummy
#define __bend_calc_m(x)\
-(1.0f - dot_rme)*(n_r.x - dot_rme*n_rme.x)*n_rme.w\
+ (1.0f - dot_l)*(n_rme.x - dot_l*n_l.x)*n_l.w - (n_l.x - dot_l*n_rme.x)*n_rme.w\
+ (1.0f - dot_ll)*(n_ll.x - dot_ll*n_l.x)*n_l.w
template<>
__device__ inline real bend_calc_1<0>(const real dot_rme, const real dot_l, const  real dot_ll, const  real4 n_r, const real4 n_rme, const real4 n_l, const  real4 n_ll){
	return __bend_calc_m(x);
}

template<>
__device__ inline real bend_calc_1<1>(const real dot_rme, const real dot_l, const  real dot_ll, const  real4 n_r, const real4 n_rme, const real4 n_l, const  real4 n_ll){
	return __bend_calc_m(y);
}
template<>
__device__ inline real bend_calc_1<2>(const real dot_rme, const real dot_l, const  real dot_ll, const  real4 n_r, const real4 n_rme, const real4 n_l, const  real4 n_ll){
	return __bend_calc_m(z);
}

__device__ inline real4 bend_calc_main(const CellPos cme, const CellPos cu, const CellPos cl, const CellPos cb, const CellPos cr,
	const CellPos cuu, const CellPos cll, const CellPos cbb, const CellPos crr){
	const real4 n_u = normal_with_len_inv(cuu, cu); //4th element == len_inv
	const real4 n_ume = normal_with_len_inv(cu, cme);
	const real4 n_b = normal_with_len_inv(cme, cb);
	const real4 n_bb = normal_with_len_inv(cb, cbb);

	const real4 n_r = normal_with_len_inv(crr, cr);
	const real4 n_rme = normal_with_len_inv(cr, cme);
	const real4 n_l = normal_with_len_inv(cme, cl);
	const real4 n_ll = normal_with_len_inv(cl, cll);

	const real dot_ume = dot3(n_u, n_ume);
	const real dot_b = dot3(n_ume, n_b);
	const real dot_bb = dot3(n_b, n_bb);

	const real dot_rme = dot3(n_r, n_rme);
	const real dot_l = dot3(n_rme, n_l);
	const real dot_ll = dot3(n_l, n_ll);



	const real ax = bend_calc_1<0>(dot_rme, dot_l, dot_ll, n_r, n_rme, n_l, n_ll)
		+ bend_calc_1<0>(dot_ume, dot_b, dot_bb, n_u, n_ume, n_b, n_bb);
	const real ay = bend_calc_1<1>(dot_rme, dot_l, dot_ll, n_r, n_rme, n_l, n_ll)
		+ bend_calc_1<1>(dot_ume, dot_b, dot_bb, n_u, n_ume, n_b, n_bb);
	const real az = bend_calc_1<2>(dot_rme, dot_l, dot_ll, n_r, n_rme, n_l, n_ll)
		+ bend_calc_1<2>(dot_ume, dot_b, dot_bb, n_u, n_ume, n_b, n_bb);
	return make_real4(ax, ay, az,0.f);
}

__global__ void memb_bend_apply(int nmemb,const CellPos*const RESTRICT cs, CellPos*const RESTRICT new_cs){
	const CellIndex i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i < nmemb){
		const CellPos cme = cs[i];
		CellPos cu = cs[get_adj_memb_idx<DIR_U>(i)];
		CellPos cl = cs[get_adj_memb_idx<DIR_L>(i)];
		CellPos cb = cs[get_adj_memb_idx<DIR_B>(i)];
		CellPos cr = cs[get_adj_memb_idx<DIR_R>(i)];

		CellPos cuu = cs[get_adj_memb_idx<DIR_U>(get_adj_memb_idx<DIR_U>(i))];
		CellPos cll = cs[get_adj_memb_idx<DIR_L>(get_adj_memb_idx<DIR_L>(i))];
		CellPos cbb = cs[get_adj_memb_idx<DIR_B>(get_adj_memb_idx<DIR_B>(i))];
		CellPos crr = cs[get_adj_memb_idx<DIR_R>(get_adj_memb_idx<DIR_R>(i))];

		const real4 straight_a = bend_calc_main(cme, cu, cl, cb, cr, cuu, cll, cbb, crr);

		cu = cs[get_adj_memb_idx<DIR_L>(get_adj_memb_idx<DIR_U>(i))];
		cl = cs[get_adj_memb_idx<DIR_B>(get_adj_memb_idx<DIR_L>(i))];
		cb = cs[get_adj_memb_idx<DIR_R>(get_adj_memb_idx<DIR_B>(i))];
		cr = cs[get_adj_memb_idx<DIR_U>(get_adj_memb_idx<DIR_R>(i))];



		cuu = cs[get_adj_memb_idx<DIR_L>(get_adj_memb_idx<DIR_U>(get_adj_memb_idx<DIR_L>(get_adj_memb_idx<DIR_U>(i))))];
		cll = cs[get_adj_memb_idx<DIR_B>(get_adj_memb_idx<DIR_L>(get_adj_memb_idx<DIR_B>(get_adj_memb_idx<DIR_L>(i))))];
		cbb = cs[get_adj_memb_idx<DIR_R>(get_adj_memb_idx<DIR_B>(get_adj_memb_idx<DIR_R>(get_adj_memb_idx<DIR_B>(i))))];
		crr = cs[get_adj_memb_idx<DIR_U>(get_adj_memb_idx<DIR_R>(get_adj_memb_idx<DIR_U>(get_adj_memb_idx<DIR_R>(i))))];

		const real4 diag_a = bend_calc_main(cme, cu, cl, cb, cr, cuu, cll, cbb, crr);


		new_cs[i] = cme + (DT_Cell*KBEND)*(straight_a + diag_a);
	}
}

__device__ real CI_memb_to_memb(const CellPos c1, const CellPos c2) {



	//using namespace cont;
	const real rad_sum = R_memb+R_memb;
	const real rad_sum_sq = rad_sum*rad_sum;

	real cr_dist_sq;
	const real distSq = p_cell_dist_sq(c1, c2);
	if (distSq< (cr_dist_sq = rad_sum_sq*P_MEMB*P_MEMB)) {
		const real LJ6 = POW3(cr_dist_sq) / POW3(distSq);
		return 4.0f*eps_m*LJ6*(LJ6 - 1.0f) / distSq;
	}
	else if (distSq < rad_sum_sq) {
		const real inv_distlj = rsqrtr(distSq);
		const real cr_dist = rad_sum*P_MEMB;
		return -(DER_DER_CONST) * (1.0f/ cr_dist - inv_distlj);
	}
	else {
		const real distlj = sqrtf(distSq);
		const real lambda_dist = (1.0f + P_MEMB)*rad_sum;
		const real cr_dist = rad_sum*P_MEMB;
		const real LJ6 = POW6(cr_dist) / POW6(lambda_dist - distlj);
		return -(DER_DER_CONST / rad_sum)*((1.0f - P_MEMB) / P_MEMB)
			- 4.0f * eps_m*(LJ6*(LJ6 - 1.0f)) / ((lambda_dist - distlj)*distlj);
	}
}

__device__ real CI_memb_to_memb_diag(int myidx,int pidx,const CellPos c1, const CellPos c2) {



	//using namespace cont;
	const real rad_sum = 2.0*R_memb;
	const real rad_sum_sq = 4.0*R_memb*R_memb;

	real cr_dist_sq;
	const real distSq = p_cell_dist_sq(c1, c2);
	if (distSq< (cr_dist_sq = rad_sum_sq*P_MEMB*P_MEMB)) {
		const real LJ6 = POW3(cr_dist_sq) / POW3(distSq);
		return 4.0f*eps_m*LJ6*(LJ6 - 1.0f) / distSq;
	}
	return 0.0f;
}
__device__ inline real ljmain(const CellPos c1, const CellPos c2,const real rad1,const real rad2) {

	const real distSq = p_cell_dist_sq(c1, c2);
	const real rad_sum = rad1 + rad2;
	const real LJ6 = POW6(rad_sum) / POW3(distSq);
	//LJ6 = LJ6*LJ6*LJ6;
	return 4.0f*eps_m*LJ6*(LJ6 - 1.0f) / distSq;

}

__device__ inline real diff_wall(const CellPos c1,const real rad1){
	const real distlj = 2.0f*c1.z;
	const real LJ6 = POW6(rad1) / POW6(c1.z);
	const real ljm = 4.0f*eps_m*LJ6*(LJ6 - 1.0f) / (distlj*distlj);
	return ljm*2.0f*c1.z;
}

__device__ inline real CI_other(const CellPos c1, const CellPos c2, real rad1, real rad2){
	return ljmain(c1, c2, rad1, rad2);
}

__global__ void memb_interaction_internal(int nmemb,
	const CellPos*const cs,
	CellPos*const new_cs,
	const CellConnectionData*const cis,
	const CELL_STATE*const cstate,
	const CellIndex*const dermis_index){
	const CellIndex i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i < nmemb){
		real4 vel = { 0.0f };
		const CellPos cme = cs[i];
		for (int j = 0; j < 4; j++){
			const CellPos opcs = cs[cis[i].connect_index[j]];
			vel = vel + CI_memb_to_memb(cme, opcs)*p_diff4(cme, opcs);
		}

		for (int j = 4; j < 8; j++){
			const CellPos opcs = cs[cis[i].connect_index[j]];
			vel = vel + CI_memb_to_memb_diag(i, cis[i].connect_index[j], cme, opcs)*p_diff4(cme, opcs);
		}

		for (int j = 8; j < cis[i].connect_num; j++){
			const CellIndex opidx = cis[i].connect_index[j];
			if (i != dermis_index[opidx]){
				const CellPos opcs = cs[opidx];
				vel = vel + CI_other(cme, opcs, R_memb, get_radius(cstate[opidx]))*p_diff4(cme, opcs);
			}
		}

		vel.z += diff_wall(cme, R_memb);
		new_cs[i] = new_cs[i] + DT_Cell*vel; //need add
	}
}

void memb_interacion(CellManager* cm){
	memb_bend_apply << <cm->nmemb_host / 256 + 1, 256 >> >(cm->nmemb_host, cm->current_pos_host(), cm->next_pos_host());
	cudaDeviceSynchronize();
	
	memb_interaction_internal << <cm->nmemb_host / 256 + 1, 256 >> >(
		cm->nmemb_host,
		cm->current_pos_host(),
		cm->next_pos_host(),
		cm->connection_data,
		cm->state,
		cm->dermis_index);
		
	cudaDeviceSynchronize();
}

__device__ inline bool is_der_near(const CellPos c1, const CellPos c2) {
	return sqrtf(p_cell_dist_sq(c1, c2)) + delta_R < 2.0f*R_der;
}
__device__ inline real ljmain_der_near(const CellPos c1, const CellPos c2) {
	/*
	‹——£‚Ì2æ‚¾‚¯‚Å‚Í•s‰Â
	*/
	const real dist = sqrtf(p_cell_dist_sq(c1, c2));
	const real dist2 = dist + delta_R;
	const real rad_sum = 2.0f*R_der;
	const real LJ6 = POW6(rad_sum) / POW6(dist2);
	//LJ6 = LJ6*LJ6*LJ6;
	return 4.0f*eps_m*LJ6*(LJ6 - 1.0f) / (dist*dist2) + para_ljp2;
}

__device__ inline bool is_near(const CellPos c1, const CellPos c2,const real rad1,const real rad2) {
	//double rad_sum = c1->radius + c2->radius;
	return p_cell_dist_sq(c1, c2) < (rad1 + rad2)*(rad1 + rad2);
}
__device__ inline real ljmain_der_far(const CellPos c1, const CellPos c2) {
	const real dist = sqrtf(p_cell_dist_sq(c1, c2));
	const real dist2 = dist + delta_R;
	const real rad_sum = 2.0f*R_der;
	const real LJ6 = POW6(rad_sum) / POW6(dist2);
	return 4.0*eps_m*LJ6*(LJ6 - 1.0f) / (dist*dist) + para_ljp2;
}

__device__ real CI_der_and_der(const CellPos c1, const CellPos c2){
	if (is_der_near(c1, c2)) {
		return ljmain_der_near(c1, c2);
	}
	else if (is_near(c1, c2, R_der, R_der)) {
		return para_ljp2;
	}
	else {
		return ljmain_der_far(c1, c2);
	}
}

__device__ real4 der_interact_impl(const CellPos cme, const CellConnectionData*const conn_ptr, const CellPos*const RESTRICT cs, const CELL_STATE*const RESTRICT cstate){
	real4 vel = { 0.0f };

	//cell_pos_set cme = cs[index];
	for (int j = 0; j < conn_ptr->connect_num; j++){
		const CellIndex opidx = conn_ptr->connect_index[j];
		const CellPos opcs = cs[opidx];
		const CELL_STATE opstate = cstate[opidx];

		const real coef = opstate == DER ? CI_der_and_der(cme, opcs) : CI_other(cme, opcs, R_der, get_radius(opstate));
		vel = vel + coef*p_diff4(cme, opcs);
	}
	vel.z += diff_wall(cme, R_der);
	return cme + DT_Cell*vel;
}

__device__ inline real adhesion(const CellPos c1,const CellPos c2,const real rad1,const real rad2, const real spring_const) {
	//using namespace cont;
	const real distlj = sqrtf(p_cell_dist_sq(c1, c2));
	const real rad_sum = rad1 + rad2;
	//double LJ2_m1;

	const real LJ2_m1 = (distlj / rad_sum) - 1.0f;
	return (LJ2_m1 + 1 > LJ_THRESH ?
		0.0f :
		-(spring_const / distlj)
		* LJ2_m1*(1.0f - LJ2_m1*LJ2_m1 / ((LJ_THRESH - 1.0f)*(LJ_THRESH - 1.0f))));
}

__device__ inline real CI_fix_mu_and_fix_mu(const CellPos c1, const CellPos c2) {
	//pair‚©‚Ç‚¤‚©‚ÍŠO‚Å”»’è‚µ‚æ‚¤
	if (is_near(c1, c2, R_max,R_max)) {
		return ljmain(c1, c2, R_max, R_max);
	}
	else {
		return adhesion(c1, c2, R_max, R_max, K_DESMOSOME);
	}
}

//always fix -> memb
__device__ real CI_fix_and_memb(const CellPos c1, const CellPos c2) {
	if (is_near(c1, c2, R_max, R_memb)) {
		return ljmain(c1, c2, R_max, R_memb);
	}
	else {
		const real distlj = sqrtf(p_cell_dist_sq(c1, c2));
		const real LJ2 = distlj / (R_max+ R_memb);
		return -(Kspring / distlj)*(LJ2 - 1.0f);
	}
}

__device__ real CI_al_air_de_fix_mu_and_al_air_de(const CellPos c1, const CellPos c2, real agek1, real agek2) {
	if (is_near(c1, c2, R_max, R_max)) {
		return ljmain(c1, c2, R_max, R_max);
	}
	else {
		const real spf = agek1 > THRESH_SP && agek2 > THRESH_SP? K_TOTAL:K_DESMOSOME;
		return adhesion(c1, c2, R_max, R_max, spf);
	}
}

__device__ real4 fix_interact_impl(const CellPos cme, const CellConnectionData*const conn_ptr, const CellPos*const RESTRICT cs, CellPos*const RESTRICT new_cs,//rewrite cs
	const CELL_STATE*const RESTRICT cstate,const real*const agek_arr,const real my_agek,const CellIndex pair_idx,const CellIndex dermis_idx){
	real4 vel = { 0.0f };
	//cell_pos_set cme = cs[index];
	//const int pair_idx = pair_arr[index];
	//const int dermis_idx = dermis_arr[index];
	for (int j = 0; j < conn_ptr->connect_num; j++){
		const int opidx = conn_ptr->connect_index[j];
		if (pair_idx == opidx)continue;


		const CellPos opcs = cs[opidx];
		const CELL_STATE opstate = cstate[opidx];
		const real4 vs = p_diff4(cme, opcs);
		real coef = 0.0f;
		switch (opstate){
		case FIX:case MUSUME:
			coef = CI_fix_mu_and_fix_mu(cme, opcs);
			break;
		case MEMB:
			if (dermis_idx == opidx){
				coef = CI_fix_and_memb(cme, opcs);
				atomicAdd(&new_cs[opidx].x, -DT_Cell*coef*vs.x);
				atomicAdd(&new_cs[opidx].y, -DT_Cell*coef*vs.y);
				atomicAdd(&new_cs[opidx].z, -DT_Cell*coef*vs.z);
			}
			else{
				coef = CI_other(cme, opcs, R_max, R_memb);
			}
			break;
		case ALIVE:case AIR:case DEAD:
			coef = CI_al_air_de_fix_mu_and_al_air_de(cme, opcs, my_agek, agek_arr[opidx]);
			break;
		case DER:
			coef = CI_other(cme, opcs, R_max, R_max);
			break;
		default:
			break;
		}
		vel = vel + coef*vs;
	}
	return cme + DT_Cell*vel;
}
//always mu->memb
__device__ real CI_mu_and_memb(const CellPos c1, const CellPos c2,const real spring_force_to_memb) {
	//spring_force_to_memb‚ªÝ’è‚³‚ê‚Ä‚¢‚é•K—v‚ª‚ ‚é
	if (is_near(c1, c2, R_max,R_memb)) {
		return ljmain(c1, c2, R_max, R_memb);
	}
	else {
		return adhesion(c1, c2, R_max, R_memb, spring_force_to_memb);
	}
}

__device__ real4 musume_interact_impl(const CellPos cme, const CellConnectionData*const conn_ptr, const CellPos*const RESTRICT cs, CellPos*const RESTRICT new_cs, //rewrite cs
	const CELL_STATE*const RESTRICT cstate, const real*const agek_arr, const int my_rest_div_time,const real my_agek, const CellIndex pair_idx, const CellIndex dermis_idx){
	real4 vel = { 0.0f };
	//cell_pos_set cme = cs[index];
	//const int pair_idx = pair_arr[index];
	//const int dermis_idx = dermis_arr[index];


	real spring_force_to_memb = 0.0f;
	if (pair_idx >= 0){
		if (cstate[pair_idx] == FIX) {
			spring_force_to_memb = Kspring;
		} else if(my_rest_div_time > 0){
			spring_force_to_memb = Kspring_d;
		}
	}
	else if (my_rest_div_time > 0) {
		spring_force_to_memb = Kspring_d;
	}

	for (int j = 0; j < conn_ptr->connect_num; j++){
		const int opidx = conn_ptr->connect_index[j];
		if (pair_idx == opidx)continue;


		const CellPos opcs = cs[opidx];
		const CELL_STATE opstate = cstate[opidx];
		const real4 vs = p_diff4(cme, opcs);
		real coef = 0.0f;
		switch (opstate){
		case FIX:case MUSUME:
			coef = CI_fix_mu_and_fix_mu(cme, opcs);
			break;
		case MEMB:
			if (dermis_idx == opidx){
				coef = CI_mu_and_memb(cme, opcs,spring_force_to_memb);
				atomicAdd(&new_cs[opidx].x, -DT_Cell*coef*vs.x);
				atomicAdd(&new_cs[opidx].y, -DT_Cell*coef*vs.y);
				atomicAdd(&new_cs[opidx].z, -DT_Cell*coef*vs.z);
			}
			else{
				coef = CI_other(cme, opcs, R_max, R_memb);
			}
			break;
		case ALIVE:case AIR:case DEAD:
			coef = CI_al_air_de_fix_mu_and_al_air_de(cme, opcs, my_agek, agek_arr[opidx]);
			break;
		case DER:
			coef = CI_other(cme, opcs, R_max, R_max);
			break;
		default:
			break;
		}
		if (fabs(coef) > 100)printf("err in musume_int\n");
		vel = vel + coef*vs;
	}
	return cme + DT_Cell*vel;
}

__device__ real4 al_air_de_interact_impl(const CellPos cme, const CellConnectionData*const conn_ptr,const CellPos*const RESTRICT cs, //rewrite cs
	const CELL_STATE*const RESTRICT cstate, const real*const agek_arr,  const real my_agek ){
	real4 vel = { 0.0f };
	//cell_pos_set cme = cs[index];
	for (int j = 0; j < conn_ptr->connect_num; j++){
		const int opidx = conn_ptr->connect_index[j];
		const CellPos opcs = cs[opidx];
		const unsigned int opstate = cstate[opidx];
		//const real3 vs = p_vec_sub(cme, opcs);
		real coef = 0.0f;
		switch (opstate){
		case ALIVE:case AIR:case DEAD:case FIX:case MUSUME:
			coef = CI_al_air_de_fix_mu_and_al_air_de(cme, opcs, my_agek, agek_arr[opidx]);
			break;
		default:
			coef = CI_other(cme, opcs, R_max, get_radius(opstate));
			break;
		}
		if (fabs(coef) > 100)printf("err in alair_int\n");
		vel = vel + coef*p_diff4(cme, opcs);
	}

	return cme + DT_Cell*vel;
}

__device__ real CI_pair(const CellPos c1,const CellPos c2,const real spr_nat_len) {
	const real inv_dist = rsqrtr(p_cell_dist_sq(c1, c2));
	return -Kspring_division*(1.0f - spr_nat_len*inv_dist);
	//return -force / dist;
}

__device__ void pair_interact(const CellIndex index,const CellPos*const RESTRICT cs,CellPos*const RESTRICT new_cs,
	const CellIndex pair_index,const real spr_nat_len){
	const CellPos cme = cs[index];
	const CellPos cpair = cs[pair_index];
	const real4 vs = p_diff4(cme, cpair);
	const real coef = CI_pair(cme, cpair, spr_nat_len);
	new_cs[index] = new_cs[index] + (DT_Cell*coef)*vs;
	new_cs[pair_index] = new_cs[pair_index] - (DT_Cell*coef)*vs;
}

__global__ void non_memb_interact(
	int nmemb, 
	int ncell, 
	const CellPos*const RESTRICT cs,
	CellPos*const RESTRICT new_cs,
	const CellConnectionData*const RESTRICT cis,
	const CELL_STATE*const RESTRICT cstate,
	const CellIndex*const RESTRICT dermis_index_arr,
	const CellIndex*const RESTRICT pair_arr,
	const real*const RESTRICT agek_arr,
	const int*const RESTRICT div_time_arr,
	const real*const RESTRICT spr_nat_len_arr){

	const int i = blockIdx.x*blockDim.x + threadIdx.x + nmemb;
	//__shared__ cell_pos_set cme_c[128];
	if (i < ncell){
		//cme_c[threadIdx.x] = cs[i];
		//__syncthreads();
		const CellPos cme = cs[i];
		const unsigned int state = cstate[i];
		const CellConnectionData* conn_ptr = &(cis[i]);
		const CellIndex pair_idx = pair_arr[i];
		real4 result = { 0.0f };
		switch (state){
		case DER:
			result = der_interact_impl(cme, conn_ptr, cs, cstate);
			break;
		case FIX:
			result = fix_interact_impl(cme, conn_ptr, cs, new_cs, cstate, agek_arr, agek_arr[i], pair_idx, dermis_index_arr[i]);
			break;
		case MUSUME:
			result = musume_interact_impl(cme, conn_ptr, cs, new_cs, cstate, agek_arr, agek_arr[i], div_time_arr[i], pair_idx, dermis_index_arr[i]);
			break;
		case ALIVE:case AIR:case DEAD:
			result = al_air_de_interact_impl(cme, conn_ptr, cs, cstate,agek_arr,agek_arr[i]);
			break;
		default:
			break;
		}
		__syncthreads();
		new_cs[i] = result;
		__syncthreads();

		if (pair_idx != -1 && i>pair_idx){
			pair_interact( i,cs, new_cs,pair_idx, spr_nat_len_arr[i]);
		}
	}
}

__global__ void cell_pos_periodic_fix_impl(CellPos*const pos){
	int index = blockIdx.x*blockDim.x + threadIdx.x;
	const CellPos cme = pos[index];
	if (cme.x > LX) {
		pos[index].x -= LX;
	}
	else if (cme.x < 0.0f) {
		pos[index].x += LX;
	}

	if (cme.y > LY) {
		pos[index].y -= LY;
	}
	else if (cme.y < 0.0f) {
		pos[index].y += LY;
	}
}
void cell_pos_periodic_fix(CellManager* cm){
	cell_pos_periodic_fix_impl << <cm->ncell_host / 512 + 1, 512 >> >(cm->current_pos_host());
	cudaDeviceSynchronize();
}

void cell_interact(CellManager* cm){
	memb_interacion(cm);
	
		non_memb_interact << <(cm->ncell_host - cm->nmemb_host) / 128 + 1, 128 >> >(
		cm->nmemb_host,
		cm->ncell_host,
		cm->current_pos_host(),
		cm->next_pos_host(),
		cm->connection_data,
		cm->state,
		cm->dermis_index,
		cm->pair_index,
		cm->agek,
		cm->rest_div_times,
		cm->spr_nat_len
		);
		
		
	cudaDeviceSynchronize();
	
	
}