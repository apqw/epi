#include "global_cuda.h"
#include <device_functions.h>
#include <stdio.h>
#include "utils_cuda.cuh"
#include "cell_interaction_cuda.h"
#define POW6(x) ((x)*(x)*(x)*(x)*(x)*(x))
#define POW3(x) ((x)*(x)*(x))
/** LJƒ|ƒeƒ“ƒVƒƒƒ‹‚ÌŒW” */
__constant__ static const float eps_m = 0.01f;

/*
MEMB_calc
*/
__constant__ static const float P_MEMB = 1.0f / COMPRESS_FACTOR;

/** L‚Ñ’e«ŒW” */
__constant__ static const float DER_DER_CONST = 0.05f;
__constant__ static const float K_TOTAL = 3.0f;
__constant__ static const float K_DESMOSOME_RATIO = 0.01f;
__constant__ static const float K_DESMOSOME = 3.0f*0.01f;
__constant__ static const float Kspring = 25.0f;
__constant__ static const float Kspring_d = 5.0f;

/** ‹È‚°’e«ŒW” */
__constant__ static const float KBEND = 0.25f;//0.5->0.25

/*
DER_calc
*/
__constant__ static const float delta_R = 0.4f*R_der;

__constant__ static const float para_ljp2 = 0.005f;
__constant__ static const float Kspring_division = 5.0f;
template<unsigned int d>
__device__ inline int get_adj_memb_idx_d(int my_memb_idx){
    return 0.0f;
}//dummy

template<>
__device__ inline  int get_adj_memb_idx_d<(unsigned int)U>(int my_memb_idx){
    int kk = my_memb_idx / NMX;
    return ((kk + 1) % NMY)*NMX + my_memb_idx%NMX;
}

template<>
__device__ inline int get_adj_memb_idx_d<(unsigned int)L>(int my_memb_idx){
    int jj = my_memb_idx%NMX;
    return ((int)(my_memb_idx / NMX))*NMX + (jj - 1+NMX) % NMX;
}

template<>
__device__ inline int get_adj_memb_idx_d<(unsigned int)B>(int my_memb_idx){
    int kk = my_memb_idx / NMX;
    return ((kk - 1+NMY) % NMY)*NMX + my_memb_idx%NMX;
}

template<>
__device__ inline int get_adj_memb_idx_d<(unsigned int)R>(int my_memb_idx){
    int jj = my_memb_idx%NMX;
    return ((int)(my_memb_idx / NMX))*NMX + (jj + 1) % NMX;
}

__device__ inline float3 p_vec_sub(cell_pos_set a, cell_pos_set b){
	return make_float3(p_diff_x(a.x, b.x), p_diff_y(a.y, b.y), a.z - b.z);
}

__device__ inline float4 normal_with_len_inv(cell_pos_set a, cell_pos_set b){
    float diffx = p_diff_x(a.x, b.x);
    float diffy = p_diff_y(a.y, b.y);
    float diffz = a.z - b.z;
    float len_inv= rsqrtf(diffx*diffx + diffy*diffy + diffz*diffz);
    return make_float4(diffx *len_inv, diffy *len_inv, diffz *len_inv, len_inv);
}

__device__ inline float dot(cell_pos_set a, cell_pos_set b){
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

template<unsigned int AXIS>
__device__ inline float bend_calc_1(float dot_rme, float dot_l, float dot_ll, float4 n_r, float4 n_rme, float4 n_l, float4 n_ll){
	return 0.0f;
}//dummy

#define __bend_calc_m(x)\
-(1.0f - dot_rme)*(n_r.x - dot_rme*n_rme.x)*n_rme.w\
+ (1.0f - dot_l)*(n_rme.x - dot_l*n_l.x)*n_l.w - (n_l.x - dot_l*n_rme.x)*n_rme.w\
+ (1.0f - dot_ll)*(n_ll.x - dot_ll*n_l.x)*n_l.w
template<>
__device__ inline float bend_calc_1<0>(float dot_rme, float dot_l, float dot_ll, float4 n_r, float4 n_rme, float4 n_l, float4 n_ll){
	return __bend_calc_m(x);
}

template<>
__device__ inline float bend_calc_1<1>(float dot_rme, float dot_l, float dot_ll, float4 n_r, float4 n_rme, float4 n_l, float4 n_ll){
	return __bend_calc_m(y);
}
template<>
__device__ inline float bend_calc_1<2>(float dot_rme, float dot_l, float dot_ll, float4 n_r, float4 n_rme, float4 n_l, float4 n_ll){
	return __bend_calc_m(z);
}

__device__ inline float3 bend_calc_main(cell_pos_set cme, cell_pos_set cu, cell_pos_set cl, cell_pos_set cb, cell_pos_set cr,
	cell_pos_set cuu, cell_pos_set cll, cell_pos_set cbb, cell_pos_set crr){
	const float4 n_u = normal_with_len_inv(cuu, cu); //4th element == len_inv
	const float4 n_ume = normal_with_len_inv(cu, cme);
	const float4 n_b = normal_with_len_inv(cme, cb);
	const float4 n_bb = normal_with_len_inv(cb, cbb);

	const float4 n_r = normal_with_len_inv(crr, cr);
	const float4 n_rme = normal_with_len_inv(cr, cme);
	const float4 n_l = normal_with_len_inv(cme, cl);
	const float4 n_ll = normal_with_len_inv(cl, cll);

	const float dot_ume = dot(n_u, n_ume);
	const float dot_b = dot(n_ume, n_b);
	const float dot_bb = dot(n_b, n_bb);

	const float dot_rme = dot(n_r, n_rme);
	const float dot_l = dot(n_rme, n_l);
	const float dot_ll = dot(n_l, n_ll);



	const float ax = bend_calc_1<0>(dot_rme, dot_l, dot_ll, n_r, n_rme, n_l, n_ll)
		+ bend_calc_1<0>(dot_ume, dot_b, dot_bb, n_u, n_ume, n_b, n_bb);
	const float ay = bend_calc_1<1>(dot_rme, dot_l, dot_ll, n_r, n_rme, n_l, n_ll)
		+ bend_calc_1<1>(dot_ume, dot_b, dot_bb, n_u, n_ume, n_b, n_bb);
	const float az = bend_calc_1<2>(dot_rme, dot_l, dot_ll, n_r, n_rme, n_l, n_ll)
		+ bend_calc_1<2>(dot_ume, dot_b, dot_bb, n_u, n_ume, n_b, n_bb);
	return make_float3(ax, ay, az);
}

__global__ void _memb_bend_calc(int nmemb, cell_pos_set* cs,cell_pos_set* new_cs){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i < nmemb){
        
    }
}



__device__ float CI_memb_to_memb(cell_pos_set c1, cell_pos_set c2,float rad1,float rad2) {



	using namespace cont;
	const float rad_sum = rad1+rad2;
	const float rad_sum_sq = rad_sum*rad_sum;

	float cr_dist_sq;
	const float distSq = p_cell_dist_sq_d(c1, c2);
	if (distSq< (cr_dist_sq = rad_sum_sq*P_MEMB*P_MEMB)) {
		//assert(fabs(ljmain(me, oppo)) < 1000);
		const float LJ6 = POW3(cr_dist_sq) / POW3(distSq);
		return 4.0f*eps_m*LJ6*(LJ6 - 1.0f) / distSq;
	}
	else if (distSq < rad_sum_sq) {
		const float distlj = sqrtf(distSq);
		const float cr_dist = rad_sum*P_MEMB;
		return -(DER_DER_CONST / distlj) * (distlj / cr_dist - 1.0f);
	}
	else {
		const float distlj = sqrtf(distSq);
		const float lambda_dist = (1.0f + P_MEMB)*rad_sum;
		const float cr_dist = rad_sum*P_MEMB;
		const float LJ6 = POW6(cr_dist) / POW6(lambda_dist - distlj);
		return -(DER_DER_CONST / rad_sum)*((1.0f - P_MEMB) / P_MEMB)
			- 4.0f * eps_m*(LJ6*(LJ6 - 1.0f)) / ((lambda_dist - distlj)*distlj);
	}
}

__device__ float CI_memb_to_memb_diag(cell_pos_set c1, cell_pos_set c2, float rad1, float rad2) {



	using namespace cont;
	const float rad_sum = rad1 + rad2;
	const float rad_sum_sq = rad_sum*rad_sum;

	float cr_dist_sq;
	const float distSq = p_cell_dist_sq_d(c1, c2);
	if (distSq< (cr_dist_sq = rad_sum_sq*P_MEMB*P_MEMB)) {
		//assert(fabs(ljmain(me, oppo)) < 1000);
		const float LJ6 = POW3(cr_dist_sq) / POW3(distSq);
		return 4.0f*eps_m*LJ6*(LJ6 - 1.0f) / distSq;
	}
	return 0.0f;
}

__device__ inline float ljmain(cell_pos_set c1, cell_pos_set c2, float rad1, float rad2) {

	const float distSq = p_cell_dist_sq_d(c1, c2);
	const float rad_sum = rad1+rad2;
	const float LJ6 = POW6(rad_sum) / POW3(distSq);
	//LJ6 = LJ6*LJ6*LJ6;
	return 4.0f*eps_m*LJ6*(LJ6 - 1.0f) / distSq;

}

__device__ inline float CI_other(cell_pos_set c1, cell_pos_set c2, float rad1, float rad2){
	return ljmain(c1, c2, rad1, rad2);
}

__device__ inline float diff_wall(cell_pos_set c1,float rad1){
	const float distlj = 2.0f*c1.z;
	const float LJ6 = POW6(rad1) / POW6(c1.z);
	const float ljm = 4.0f*eps_m*LJ6*(LJ6 - 1.0f) / (distlj*distlj);
	return ljm*2.0f*c1.z;
}

__global__ void memb_bend_apply_2(int nmemb, cell_pos_set* cs, cell_pos_set* new_cs, connected_index_set* cis){
	const int i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i < nmemb){
		const cell_pos_set cme = cs[i];
		cell_pos_set cu = cs[get_adj_memb_idx_d<U>(i)];
		cell_pos_set cl = cs[get_adj_memb_idx_d<L>(i)];
		cell_pos_set cb = cs[get_adj_memb_idx_d<B>(i)];
		cell_pos_set cr = cs[get_adj_memb_idx_d<R>(i)];

		cell_pos_set cuu = cs[get_adj_memb_idx_d<U>(get_adj_memb_idx_d<U>(i))];
		cell_pos_set cll = cs[get_adj_memb_idx_d<L>(get_adj_memb_idx_d<L>(i))];
		cell_pos_set cbb = cs[get_adj_memb_idx_d<B>(get_adj_memb_idx_d<B>(i))];
		cell_pos_set crr = cs[get_adj_memb_idx_d<R>(get_adj_memb_idx_d<R>(i))];

		const float3 straight_a = bend_calc_main(cme, cu, cl, cb, cr, cuu, cll, cbb, crr);

		cu = cs[get_adj_memb_idx_d<L>(get_adj_memb_idx_d<U>(i))];
		cl = cs[get_adj_memb_idx_d<B>(get_adj_memb_idx_d<L>(i))];
		cb = cs[get_adj_memb_idx_d<R>(get_adj_memb_idx_d<B>(i))];
		cr = cs[get_adj_memb_idx_d<U>(get_adj_memb_idx_d<R>(i))];



		cuu = cs[get_adj_memb_idx_d<L>(get_adj_memb_idx_d<U>(get_adj_memb_idx_d<L>(get_adj_memb_idx_d<U>(i))))];
		cll = cs[get_adj_memb_idx_d<B>(get_adj_memb_idx_d<L>(get_adj_memb_idx_d<B>(get_adj_memb_idx_d<L>(i))))];
		cbb = cs[get_adj_memb_idx_d<R>(get_adj_memb_idx_d<B>(get_adj_memb_idx_d<R>(get_adj_memb_idx_d<B>(i))))];
		crr = cs[get_adj_memb_idx_d<U>(get_adj_memb_idx_d<R>(get_adj_memb_idx_d<U>(get_adj_memb_idx_d<R>(i))))];

		const float3 diag_a = bend_calc_main(cme, cu, cl, cb, cr, cuu, cll, cbb, crr);


		new_cs[i].x = cme.x + cont::DT_Cell*KBEND*(straight_a.x + diag_a.x);
		new_cs[i].y = cme.y + cont::DT_Cell*KBEND*(straight_a.y + diag_a.y);
		new_cs[i].z = cme.z + cont::DT_Cell*KBEND*(straight_a.z + diag_a.z);
	}
}

__global__ void cell_interaction_memb_apply(int nmemb,cell_pos_set* cs, cell_pos_set* new_cs,connected_index_set* cis,unsigned int* cstate,int* dermis_index_arr) {
	const int i = blockIdx.x*blockDim.x + threadIdx.x;
	__shared__ cell_pos_set cme_c[128];
	if (i < nmemb){

		//do after bend calc
		float4 vel = { 0.0f };


		cme_c[threadIdx.x] = cs[i];

		

		

		//int opidx = -1;
		cell_pos_set opcs;
		//float4 vs;
		float coef;
#define CME cme_c[threadIdx.x]
		for (int j = 0; j < 4; j++){
			//opidx = cis[i].index[j];
			opcs = cs[cis[i].index[j]];
			//vs = p_vec_sub(CME, opcs);
			coef = CI_memb_to_memb(CME, opcs, R_memb, R_memb);
			vel.x += coef*p_diff_x(CME.x,opcs.x);
			vel.y += coef*p_diff_y(CME.y, opcs.y);
			vel.z += coef*(CME.z-opcs.z);
		}

		for (int j = 4; j < 8; j++){
			//opidx = cis[i].index[j];
			opcs = cs[cis[i].index[j]];
			//vs = p_vec_sub(CME, opcs);
			coef = CI_memb_to_memb_diag(CME, opcs, R_memb, R_memb);
			vel.x += coef*p_diff_x(CME.x, opcs.x);
			vel.y += coef*p_diff_y(CME.y, opcs.y);
			vel.z += coef*(CME.z - opcs.z);
		}
		
		
		for (int j = 8; j < cis[i].connected_num; j++){
			//opidx = cis[i].index[j];
			if (i != dermis_index_arr[cis[i].index[j]]){
				opcs = cs[cis[i].index[j]];
				//vs = p_vec_sub(CME, opcs);
				coef = CI_other(CME, opcs, R_memb, R_max);
				vel.x += coef*p_diff_x(CME.x, opcs.x);
				vel.y += coef*p_diff_y(CME.y, opcs.y);
				vel.z += coef*(CME.z - opcs.z);
			}
		}

		vel.z += diff_wall(CME, R_memb);
		
		//this is right because after bend calc
		new_cs[i].x += cont::DT_Cell*vel.x;
		new_cs[i].y += cont::DT_Cell*vel.y;
		new_cs[i].z += cont::DT_Cell*vel.z;
	}
}

__device__ inline bool is_der_near(cell_pos_set c1, cell_pos_set c2, float rad1, float rad2) {
	return sqrtf(p_cell_dist_sq_d(c1,c2)) + delta_R < rad1+rad2;
}

__device__ inline float ljmain_der_near(cell_pos_set c1, cell_pos_set c2, float rad1, float rad2) {
	/*
	‹——£‚Ì2æ‚¾‚¯‚Å‚Í•s‰Â
	*/
	const float dist = sqrtf(p_cell_dist_sq_d(c1, c2));
	const float dist2 = dist + delta_R;
	const float rad_sum = rad1 + rad2;
	const float LJ6 = POW6(rad_sum) / POW6(dist2);
	//LJ6 = LJ6*LJ6*LJ6;
	return 4.0f*eps_m*LJ6*(LJ6 - 1.0f) / (dist*dist2) + para_ljp2;
}

__device__ inline bool is_near(cell_pos_set c1, cell_pos_set c2, float rad1, float rad2) {
	//double rad_sum = c1->radius + c2->radius;
	return p_cell_dist_sq_d(c1, c2) < (rad1+rad2)*(rad1+rad2);
}

__device__ inline float ljmain_der_far(cell_pos_set c1, cell_pos_set c2, float rad1, float rad2) {
	const float dist = sqrtf(p_cell_dist_sq_d(c1, c2));
	const float dist2 = dist + delta_R;
	const float rad_sum = rad1 + rad2;
	const float LJ6 = POW6(rad_sum) / POW6(dist2);
	return 4.0*eps_m*LJ6*(LJ6 - 1.0f) / (dist*dist) + para_ljp2;
}
__device__ float CI_der_and_der(cell_pos_set c1, cell_pos_set c2, float rad1, float rad2){
	if (is_der_near(c1, c2,rad1,rad2)) {
		return ljmain_der_near(c1, c2, rad1, rad2);
	}
	else if (is_near(c1, c2, rad1, rad2)) {
		return para_ljp2;
	}
	else {
		return ljmain_der_far(c1, c2, rad1, rad2);
	}
}


__device__ float3 der_interact_impl(cell_pos_set cme, int index, cell_pos_set*__restrict cs, cell_pos_set*__restrict new_cs, connected_index_set*__restrict cis, unsigned int*__restrict cstate){
	float4 vel = { 0.0f };
	//cell_pos_set cme = cs[index];
	for (int j = 0; j < cis[index].connected_num; j++){
		int opidx = cis[index].index[j];
		cell_pos_set opcs = cs[opidx];
		unsigned int opstate = cstate[opidx];
		float3 vs = p_vec_sub(cme, opcs);
		float coef = opstate == DER ? CI_der_and_der(cme, opcs, get_radius(DER), get_radius(DER)) : CI_other(cme, opcs, get_radius(DER), get_radius(opstate));
		vel.x += coef*vs.x;
		vel.y += coef*vs.y;
		vel.z += coef*vs.z;
	}
	vel.z += diff_wall(cme, get_radius(DER));

	return make_float3(cme.x + cont::DT_Cell*vel.x, cme.y + cont::DT_Cell*vel.y, cme.z + cont::DT_Cell*vel.z);
	/*
	new_cs[index].x = cme.x + cont::DT_Cell*vel.x;
	new_cs[index].y = cme.y + cont::DT_Cell*vel.y;
	new_cs[index].z = cme.z + cont::DT_Cell*vel.z;
	*/
}

__device__ inline float adhesion(cell_pos_set c1, cell_pos_set c2, float rad1, float rad2, const float spring_const) {
	using namespace cont;
	const float distlj = sqrtf(p_cell_dist_sq_d(c1, c2));
	const float rad_sum = rad1+rad2;
	//double LJ2_m1;

	const float LJ2_m1 = (distlj / rad_sum) - 1.0f;
	return (LJ2_m1 + 1 > LJ_THRESH ?
		0.0f :
		-(spring_const / distlj)
		* LJ2_m1*(1.0f - LJ2_m1*LJ2_m1 / ((LJ_THRESH - 1.0f)*(LJ_THRESH - 1.0f))));
}

__device__ inline float CI_fix_mu_and_fix_mu(cell_pos_set c1, cell_pos_set c2, float rad1, float rad2) {
	//pair‚©‚Ç‚¤‚©‚ÍŠO‚Å”»’è‚µ‚æ‚¤
	if (is_near(c1, c2,rad1,rad2)) {
		return ljmain(c1, c2,rad1,rad2);
	}
	else {
		return adhesion(c1, c2,rad1,rad2, K_DESMOSOME);
	}
}

__device__ float CI_fix_and_memb(cell_pos_set c1, cell_pos_set c2, float rad1, float rad2) {
	if (is_near(c1, c2,rad1,rad2)) {
		return ljmain(c1, c2, rad1, rad2);
	}
	else {
		const float distlj = sqrtf(p_cell_dist_sq_d(c1, c2));
		const float LJ2 = distlj / (rad1+rad2);
		return -(Kspring / distlj)*(LJ2 - 1.0f);
	}
}

__device__ float CI_al_air_de_fix_mu_and_al_air_de(cell_pos_set c1, cell_pos_set c2, float rad1, float rad2,float agek1,float agek2) {
	if (is_near(c1, c2, rad1, rad2)) {
		return ljmain(c1, c2, rad1, rad2);
	}
	else {
		float spf = K_DESMOSOME;
		if (agek1 > cont::THRESH_SP && agek2 > cont::THRESH_SP) {
			spf = K_TOTAL;
		}
		return adhesion(c1, c2,rad1,rad2, spf);
	}
}

__device__ float3 fix_interact_impl(cell_pos_set cme, int index, cell_pos_set*__restrict cs, cell_pos_set*__restrict new_cs, connected_index_set* cis, unsigned int* cstate, int* dermis_arr, int* pair_arr, float* agek_arr){
	float4 vel = { 0.0f };
	//cell_pos_set cme = cs[index];
	const int pair_idx = pair_arr[index];
	const int dermis_idx = dermis_arr[index];
	for (int j = 0; j < cis[index].connected_num; j++){
		const int opidx = cis[index].index[j];
		if (pair_idx == opidx)continue;


		const cell_pos_set opcs = cs[opidx];
		const unsigned int opstate = cstate[opidx];
		const float3 vs = p_vec_sub(cme, opcs);
		float coef = 0.0f;
		switch (opstate){
		case FIX:case MUSUME:
			coef = CI_fix_mu_and_fix_mu(cme, opcs, get_radius(FIX), get_radius(opstate));
			break;
		case MEMB:
			if (dermis_idx == opidx){
				coef = CI_fix_and_memb(cme, opcs, get_radius(FIX), get_radius(MEMB));
				atomicAdd(&cs[opidx].x, -cont::DT_Cell*coef*vs.x);
				atomicAdd(&cs[opidx].y, -cont::DT_Cell*coef*vs.y);
				atomicAdd(&cs[opidx].z, -cont::DT_Cell*coef*vs.z);
			}
			else{
				coef = CI_other(cme, opcs, get_radius(FIX), get_radius(MEMB));
			}
			break;
		case ALIVE:case AIR:case DER:
			coef = CI_al_air_de_fix_mu_and_al_air_de(cme, opcs, get_radius(FIX), get_radius(opstate), agek_arr[index], agek_arr[opidx]);
			break;
		default:
			coef = CI_other(cme, opcs, get_radius(FIX), get_radius(opstate));
			break;
		}
		vel.x += coef*vs.x;
		vel.y += coef*vs.y;
		vel.z += coef*vs.z;
	}
	return make_float3(cme.x + cont::DT_Cell*vel.x, cme.y + cont::DT_Cell*vel.y, cme.z + cont::DT_Cell*vel.z);
	/*
	new_cs[index].x = cme.x + cont::DT_Cell*vel.x;
	new_cs[index].y = cme.y + cont::DT_Cell*vel.y;
	new_cs[index].z = cme.z + cont::DT_Cell*vel.z;
	*/
}



__device__ float CI_mu_and_memb(cell_pos_set c1, cell_pos_set c2, float rad1, float rad2, float spring_force_to_memb) {
	//spring_force_to_memb‚ªÝ’è‚³‚ê‚Ä‚¢‚é•K—v‚ª‚ ‚é
	if (is_near(c1, c2,rad1,rad2)) {
		return ljmain(c1, c2, rad1, rad2);
	}
	else {
		return adhesion(c1, c2, rad1, rad2,spring_force_to_memb);
	}
}

__device__ float3 musume_interact_impl(cell_pos_set cme, int index, cell_pos_set*__restrict cs, cell_pos_set*__restrict new_cs, connected_index_set*__restrict cis, unsigned int*__restrict cstate, int*__restrict dermis_arr, int*__restrict pair_arr, float*__restrict agek_arr, int*__restrict div_time_arr){
	float3 vel = { 0.0f };
	//cell_pos_set cme = cs[index];
	const int pair_idx = pair_arr[index];
	const int dermis_idx = dermis_arr[index];

	float spring_force_to_memb = 0;
	if (pair_idx>=0)if(cstate[pair_idx]==FIX) {
		spring_force_to_memb = Kspring;
	}
	else if (div_time_arr[index] > 0) {
		spring_force_to_memb = Kspring_d;
	}
	for (int j = 0; j < cis[index].connected_num; j++){
		const int opidx = cis[index].index[j];
		if (pair_idx == opidx)continue;


		const cell_pos_set opcs = cs[opidx];
		const unsigned int opstate = cstate[opidx];
		const float3 vs = p_vec_sub(cme, opcs);
		float coef = 0.0f;
		switch (opstate){
		case FIX:case MUSUME:
			coef = CI_fix_mu_and_fix_mu(cme, opcs, get_radius(MUSUME), get_radius(opstate));
			break;
		case MEMB:
			if (dermis_idx == opidx){
				coef = CI_mu_and_memb(cme, opcs, get_radius(MUSUME), get_radius(MEMB),spring_force_to_memb);
				atomicAdd(&cs[opidx].x, -cont::DT_Cell*coef*vs.x);
				atomicAdd(&cs[opidx].y, -cont::DT_Cell*coef*vs.y);
				atomicAdd(&cs[opidx].z, -cont::DT_Cell*coef*vs.z);
			}
			else{
				coef = CI_other(cme, opcs, get_radius(MUSUME), get_radius(MEMB));
			}
			break;
		case ALIVE:case AIR:case DER:
			coef = CI_al_air_de_fix_mu_and_al_air_de(cme, opcs, get_radius(MUSUME), get_radius(opstate), agek_arr[index], agek_arr[opidx]);
			break;
		default:
			coef = CI_other(cme, opcs, get_radius(MUSUME), get_radius(opstate));
			break;
		}
		vel.x += coef*vs.x;
		vel.y += coef*vs.y;
		vel.z += coef*vs.z;
	}

	return make_float3(cme.x + cont::DT_Cell*vel.x, cme.y + cont::DT_Cell*vel.y, cme.z + cont::DT_Cell*vel.z);
}


__device__ float3 al_air_de_interact_impl(cell_pos_set cme, int index, cell_pos_set*__restrict cs, cell_pos_set*__restrict new_cs, connected_index_set*__restrict cis, unsigned int*__restrict cstate, float*__restrict agek_arr){
	float3 vel = { 0.0f };
	//cell_pos_set cme = cs[index];
	for (int j = 0; j < cis[index].connected_num; j++){
		const int opidx = cis[index].index[j];
		const cell_pos_set opcs = cs[opidx];
		const unsigned int opstate = cstate[opidx];
		const float3 vs = p_vec_sub(cme, opcs);
		float coef = 0.0f;
		switch (opstate){
		case ALIVE:case AIR:case DEAD:case FIX:case MUSUME:
			coef = CI_al_air_de_fix_mu_and_al_air_de(cme, opcs, R_max, get_radius(opstate), agek_arr[index], agek_arr[opidx]);
			break;
		default:
			coef = CI_other(cme, opcs, get_radius(MUSUME), get_radius(opstate));
			break;
		}
		vel.x += coef*vs.x;
		vel.y += coef*vs.y;
		vel.z += coef*vs.z;
	}

	return make_float3(cme.x + cont::DT_Cell*vel.x, cme.y + cont::DT_Cell*vel.y, cme.z + cont::DT_Cell*vel.z);
}

__device__ float CI_pair(cell_pos_set c1, cell_pos_set c2,float spr_nat_len) {
	const float dist = sqrtf(p_cell_dist_sq_d(c1, c2));
	const float force = Kspring_division*(dist - spr_nat_len);
	return -force / dist;
}

__device__ void pair_interact_impl(cell_pos_set cme, int index, int pair_index, cell_pos_set*__restrict cs, cell_pos_set*__restrict new_cs, float spr_nat_len){
	//cell_pos_set cme = cs[index];
	const cell_pos_set cpair = cs[pair_index];
	float3 vs = p_vec_sub(cme, cpair);
	float coef = CI_pair(cme, cpair, spr_nat_len);
	new_cs[index].x += cont::DT_Cell*coef*vs.x;
	new_cs[index].y += cont::DT_Cell*coef*vs.y;
	new_cs[index].z += cont::DT_Cell*coef*vs.z;

	new_cs[pair_index].x -=  cont::DT_Cell*coef*vs.x;
	new_cs[pair_index].y -=  cont::DT_Cell*coef*vs.y;
	new_cs[pair_index].z -=  cont::DT_Cell*coef*vs.z;
}
__global__ void non_memb_interact_impl(int nmemb,int ncell, cell_pos_set* cs, cell_pos_set* new_cs, connected_index_set* cis, unsigned int* cstate, int* dermis_index_arr,int* pair_arr,float* agek_arr,int* div_time_arr,float* spr_nat_len_arr){
	const int i = blockIdx.x*blockDim.x + threadIdx.x+nmemb;
	__shared__ cell_pos_set cme_c[128];
	if (i < ncell){
		cme_c[threadIdx.x] = cs[i];
		__syncthreads();
		const unsigned int state = cstate[i];
		float3 result = { 0.0f };
		switch (state){
		case DER:
			result=der_interact_impl(cme_c[threadIdx.x],i, cs, new_cs, cis, cstate);
			break;
		case FIX:
			result = fix_interact_impl(cme_c[threadIdx.x], i, cs, new_cs, cis, cstate, dermis_index_arr, pair_arr, agek_arr);
			break;
		case MUSUME:
			result = musume_interact_impl(cme_c[threadIdx.x], i, cs, new_cs, cis, cstate, dermis_index_arr, pair_arr, agek_arr, div_time_arr);
			break;
		case ALIVE:case AIR:case DEAD:
			result = al_air_de_interact_impl(cme_c[threadIdx.x], i, cs, new_cs, cis, cstate, agek_arr);
			break;
		default:
			break;
		}
		__syncthreads();
		new_cs[i].x = result.x;
		new_cs[i].y = result.y;
		new_cs[i].z = result.z;
		__syncthreads();

		if (pair_arr[i] != -1 && i>pair_arr[i]){
			pair_interact_impl(cme_c[threadIdx.x], i, pair_arr[i], cs, new_cs, spr_nat_len_arr[i]);
		}
	}
}
__global__ void cell_pos_periodic_fix(cell_pos_set* next){
	int index = blockIdx.x*blockDim.x + threadIdx.x;
	cell_pos_set cme = next[index];
	if (cme.x > LX) {
		next[index].x -= LX;
	}
	else if (cme.x < 0.0f) {
		next[index].x += LX;
	}

	if (cme.y > LY) {
		next[index].y -= LY;
	}
	else if (cme.y < 0.0f) {
		next[index].y += LY;
	}
}
void memb_interact(DeviceData* d){
	
	
}

void interact(DeviceData*d){
	memb_bend_apply_2 << <300 * 150 / 128 + 1, 128 >> >(d->nmemb, d->c_pos_d[d->current], d->c_pos_d[1 - d->current], d->c_connected_index_d);
	cudaDeviceSynchronize();
	cell_interaction_memb_apply << <300 * 150 / 128 + 1, 128 >> >(d->nmemb, d->c_pos_d[d->current], d->c_pos_d[1 - d->current], d->c_connected_index_d, (unsigned int*)d->c_state_d, d->c_dermis_index_d);
	non_memb_interact_impl << <(d->ncell - d->nmemb) /128+ 1, 128 >> >(d->nmemb, d->ncell, d->c_pos_d[d->current], d->c_pos_d[1 - d->current], d->c_connected_index_d, (unsigned int*)d->c_state_d, d->c_dermis_index_d, d->c_pair_index_d, d->c_agek_d, d->c_rest_div_times_d, d->c_spr_nat_len_d);
	cudaDeviceSynchronize();
	cell_pos_periodic_fix << <d->ncell / 512 + 1, 512 >> >(d->c_pos_d[1 - d->current]);
	cudaDeviceSynchronize();
}