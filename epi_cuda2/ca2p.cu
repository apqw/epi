#include "define.h"
#include "CellManager.h"
#include "utils.h"
#include "ca2p.h"
#include "device_functions.h"
#include <cstdio>
/** カルシウムの平均化に使うイテレーション回数 */
#define Ca_ITR (int)(Ca_avg_time / DT_Ca)

/** @todo 分かりやすい命名 */
#define Kpp CDEF(0.3)
/** @todo 分かりやすい命名 */
#define dp CDEF(0.1)
__global__ void init_ca2p_map(real* air_stim, const real** diffu_map, const real* diffu, const int*const cmap1, const CellConnectionData* conn, const CELL_STATE* state, int iz_bound, const real* default_diffu_ptr){
	const int iz = threadIdx.x;
	const int ix = blockIdx.x;
	const int iy = blockIdx.y;

	const int coord = midx<NX + 1, NY + 1, NZ + 1>(ix, iy, iz);
	const CellIndex index = cmap1[coord];
	const real* tmp = default_diffu_ptr;
	real asf = 0.0f;
	if (index >= 0){
		const int conn_num = conn[index].connect_num;
		
		const CELL_STATE st = state[index];
		if (st == ALIVE || st == FIX || st == MUSUME){
			tmp = &diffu[index];
		}

		for (int i = 0; i < conn_num; i++){
			if (state[conn[index].connect_index[i]] == AIR){
				asf = 1.0f;
				break;
			}
		}
		
	}
	diffu_map[coord] = tmp;
	air_stim[coord] = asf;
}

__global__ void dead_IP3_calc(int ncell, int offset, const CellConnectionData* conn, const CELL_STATE* state, const real* IP3, real* next_IP3){
	const CellIndex index = blockIdx.x*blockDim.x + threadIdx.x + offset;
	if (index < ncell){

		if (state[index] == DEAD){

			const int conn_num = conn[index].connect_num;
			real tmp = 0.0f;
			const real my_ip3 = IP3[index];
			for (int i = 0; i < conn_num; i++){
				const CellIndex cidx = conn[index].connect_index[i];
				if (state[cidx] == ALIVE){
					tmp += IP3[cidx] - my_ip3;
				}
			}
			tmp = DT_Ca*(dp*tmp - Kpp*my_ip3);
			next_IP3[index] = my_ip3 + tmp;
		}
	}
}
/** @todo 分かりやすい命名 */
__device__ inline real calc_th(bool is_alive, real agek) {

	#define thgra (0.2f)
	#define delta_th (1.0f)
	#define thpri (1.0f)
	//using namespace cont;


	return is_alive ?
		thgra + ((thpri - thgra)*0.5f) * (1.0f + tanhf((THRESH_SP - agek) / delta_th)) :
		thpri;
}

/** @todo 分かりやすい命名 */
__device__ inline real calc_Kpa(bool is_alive, real agek) {

	#define kpa (4.0f)
	#define Kgra (kpa)
	#define delta_K (1.0f)
	#define Kpri (6.0f)

	//using namespace cont;
	return is_alive ?
		Kgra + ((Kpri - Kgra)*0.5f) * (1.0f + tanhf((THRESH_SP - agek) / delta_K)) :
		Kpri;
}

/** @todo 分かりやすい命名 */
__device__ inline real calc_IAG(bool is_alive, real agek) {
	#define delta_I (1.5f)
	#define iage_kitei (0.0f)
	return is_alive ?
		0.5f*(1.0f + tanhf((agek - THRESH_SP) / delta_I)) :
		iage_kitei;
}

/** @todo 分かりやすい命名 */
__device__ inline real fw(real diff, real w)
{
	#define wd (0.1f)
	//using namespace cont;
	return (1.0f - w) + (-1.0f + tanhf((wd - diff) / 0.1f)) / 2.0f; //epsw0 == 0.1 //done
	//-1. <-????
}


////////////////////////////////////////////////////////////////////////////////
//                                                                          
//  Ca2+の反応項に関する定義
//
////////////////////////////////////////////////////////////////////////////////

__device__ inline real ca2p_by_ext_stim(real ext_stim) {
	#define kbc (0.4f*1.2f)
	#define Hb (0.01f)
	#define Cout (1.0f)

	return kbc * ext_stim*ext_stim * Cout / (Hb + ext_stim*ext_stim);
}

__device__ inline real ca2p_into_storage(real ca2p_in_cell) {
	#define kg (0.1f)
	#define gamma (2.0f)
	return gamma * ca2p_in_cell / (kg + ca2p_in_cell);
}

__device__ inline real ER_domain1_active(real ip3) {
#define mu0 (0.567f)
	#define mu1 (0.1f)
	#define kmu (0.05f)
	return (mu0 + mu1 * ip3 / (ip3 + kmu));
}

__device__ inline real ER_domain2_active(real ca2p) {
	#define para_b (0.11f)
	#define para_bb (0.89f)
	#define para_k1 (0.7f)
	return (para_b + para_bb * ca2p / (para_k1 + ca2p));
}

/**
*  Ca2+の反応項
*  @param [in] u Ca2+濃度
*  @param [in] ER_domain3_deactive 貯蔵庫から細胞質内への放出に対する不活性効果
*  @param [in] p IP3濃度
*  @param [in] B 細胞外刺激物質濃度
*/
__device__ inline real ca2p_reaction_factor(real u, real ER_domain3_deactive, real p, real B)
{
	#define k_flux (8.1f)

	#define beta_zero (0.02f)
	#define CA_OUT (1.0f)
	#define leak_from_storage (CA_OUT*beta_zero)

	const real ca2p_from_storage = k_flux * ER_domain1_active(p) * ER_domain2_active(u) * ER_domain3_deactive;

	return
		ca2p_from_storage
		- ca2p_into_storage(u)
		+ leak_from_storage
		+ ca2p_by_ext_stim(B);



}

/** @todo 分かりやすい命名 */
__device__ inline real IP3_default_diff(real _Kpa, real a_avg, real current_IP3) {
	#define H0 (0.5f)

	//using namespace cont;
	return (_Kpa*a_avg / (H0 + a_avg) - Kpp*current_IP3);
}

/** @todo 分かりやすい命名 */
__device__ inline real ex_inert_diff(real ca2p, real current_ex_inert, real _th) {
	#define para_k2 (0.7f)
	//using namespace cont;
	return ((para_k2*para_k2 / (para_k2*para_k2 + ca2p*ca2p) - current_ex_inert) / _th);
}

__global__ void supra_calc(int ncell, int offset, CellConnectionData* conn,const CellPos* pos, const CELL_STATE* state, const real* agek,const real* IP3, real* next_IP3
	, const real* ca2p, real* next_ca2p, real* ca2p_avg, real* ex_inert, real* diffu,const real* ext_stim,const real* ATP){
#define ca2p_du (0.01f)
	const CellIndex index = blockIdx.x*blockDim.x + threadIdx.x + offset;
	if (index < ncell){
		const CELL_STATE st = state[index];
		if (st == ALIVE || st == FIX || st == MUSUME){
			diffu[index] = 0.0f;
			const CellPos cs = pos[index];
			const int3 lat = {
				(int)rintf(cs.x*inv_dx) % NX,
				(int)rintf(cs.y*inv_dy) % NY,
				(int)rintf(cs.z*inv_dz) % NZ
			};
			const real my_agek = agek[index];
			const real my_ca2p = ca2p[index];
			const real my_ip3 = IP3[index];
			const real my_ex_inert = ex_inert[index];
			const bool is_alive = st == ALIVE;
			const real _th = calc_th(is_alive, my_agek);
			const real _Kpa = calc_Kpa(is_alive, my_agek);
			const real IAGv = calc_IAG(is_alive, my_agek);
			
			real tmp_diffu = 0.0f;
			real tmp_IP3 = 0.0f;

			const int conn_num = conn[index].connect_num;
			
			for (int i = 0; i < conn_num; i++){
				const CellIndex cidx = conn[index].connect_index[i];
				const CELL_STATE cst = state[cidx];
				//printf("%d aaa\n", i);
				if (cst == ALIVE || cst == DEAD || cst == FIX || cst == MUSUME){
					tmp_diffu += conn[index].current_gj()[i] * (ca2p[cidx] - my_ca2p);//TODO:set key
					tmp_IP3 += conn[index].current_gj()[i] * (IP3[cidx] - my_ip3);
					
					if (cst == ALIVE){
						//if(i>200)printf("%d aaa\n", i);
						conn[index].current_gj()[i] = DT_Ca*fw(fabsf(ca2p[cidx] - my_ca2p), conn[index].current_gj()[i]);
					}
					
				}
				
			}
			
			
			const real my_diffu = (diffu[index] = ca2p_reaction_factor(my_ca2p, my_ex_inert, my_ip3, grid_avg8<NX + 1, NY + 1, NZ + 1>(ext_stim, lat.x, lat.y, lat.z)) + ca2p_du*IAGv*tmp_diffu);

			next_IP3[index] = my_ip3 + DT_Ca*(IP3_default_diff(_Kpa, grid_avg8<NX + 1, NY + 1, NZ + 1>(ATP, lat.x, lat.y, lat.z), my_ip3) + dp*IAGv*tmp_IP3);
			
			const real n_ca2p = my_ca2p + DT_Ca*my_diffu;
			next_ca2p[index] = n_ca2p;
			ca2p_avg[index] += n_ca2p;
			ex_inert[index] = my_ex_inert + DT_Ca*ex_inert_diff(my_ca2p, my_ex_inert, _th);
			
		}
	}
}


/** @todo 分かりやすい命名 */
__device__ inline real fa(real diffu, real A) {
	#define STIM11 (0.002f)
	#define Kaa (0.5f)
	//using namespace cont;
	return STIM11*min0(diffu) - A*Kaa;
}

__global__  void ATP_refresh(const real* ATP,real* next_ATP, const real** diffu_map, const real* air_stim_arr,const float* cmap2, int iz_bound,int shmem_border) {
#define Da (1.0f)
#define AIR_STIM (0.1f)
	
	const int iz = threadIdx.x;
	const int ix = blockIdx.x;
	const int iy = blockIdx.y;

	extern __shared__ float cmap2_line[];
	real* atp_line = (real*)&cmap2_line[shmem_border];
	//__shared__ real ext_stim_line[NZ + 1];
	const int next_iy = (iy + 1) % NY;
	const int prev_iy = (iy - 1 + NY) % NY;


	const int next_ix = (ix + 1) % NX;
	const int prev_ix = (ix - 1 + NX) % NX;

	const int next_iz = iz == NZ - 1 ? NZ - 2 : iz + 1;
	const int prev_iz = iz == 0 ? 1 : iz - 1;
	//printf("t1");

	const int center = midx<NX + 1, NY + 1, NZ + 1>(ix, iy, iz);

	const int idx1 = midx<NX + 1, NY + 1, NZ + 1>(ix, prev_iy, iz);
	const int idx2 = midx<NX + 1, NY + 1, NZ + 1>(ix, next_iy, iz);
	const int idx3 = midx<NX + 1, NY + 1, NZ + 1>(prev_ix, iy, iz);
	const int idx4 = midx<NX + 1, NY + 1, NZ + 1>(next_ix, iy, iz);
	
	cmap2_line[iz] = cmap2[center];
	const real cvalue = (atp_line[iz] = ATP[center]);
	//printf("t2");
	if (iz == iz_bound - 1){
		cmap2_line[next_iz] = cmap2[midx<NX + 1, NY + 1, NZ + 1>(ix, iy, next_iz)];
		atp_line[next_iz] = ATP[midx<NX + 1, NY + 1, NZ + 1>(ix, iy, next_iz)];
	}

	if (iz == 0){
		cmap2_line[prev_iz] = cmap2[midx<NX + 1, NY + 1, NZ + 1>(ix, iy, prev_iz)];
		atp_line[prev_iz] = ATP[midx<NX + 1, NY + 1, NZ + 1>(ix, iy, prev_iz)];
	}
	
	//printf("t3");
	next_ATP[center] = cvalue + DT_Ca*(Da*(
		cmap2_line[prev_iz] * (atp_line[prev_iz] - cvalue)
		+ cmap2_line[next_iz] * (atp_line[next_iz] - cvalue)
		+ cmap2[idx1] * (ATP[idx1] - cvalue)
		+ cmap2[idx2] * (ATP[idx2] - cvalue)
		+ cmap2[idx3] * (ATP[idx3] - cvalue)
		+ cmap2[idx4] * (ATP[idx4] - cvalue)
		)*inv_dx*inv_dx + fa(*diffu_map[center], cvalue) + air_stim_arr[center] * AIR_STIM);
	//printf("t4");
}

__global__ void init_ca2p_avg(int ncell, int offset, const CELL_STATE* cstate, real* ca2p_avg){
	const CellIndex index = blockIdx.x*blockDim.x + threadIdx.x + offset;
	if (index < ncell){
		const CELL_STATE st = cstate[index];
		if (st == ALIVE || st == FIX || st == MUSUME) {
			ca2p_avg[index] = 0.0f;
		}
	}
}

/**
*  細胞にカルシウム濃度の平均値をセットする。
*/
__global__ void set_cell_ca2p(int ncell, int offset, real* ca2p_avg, const CELL_STATE* state, real* ca2p,real* next_ca2p){
	const CellIndex index = blockIdx.x*blockDim.x + threadIdx.x + offset;
	if (index < ncell){
		CELL_STATE st = state[index];
		if (st == ALIVE || st == FIX || st == MUSUME){
			ca2p_avg[index] /= Ca_ITR;
		}
		else{
			ca2p_avg[index] = 0.0f;
			ca2p[index] = 0.0f;
			next_ca2p[index] = 0.0f;
		}
	}
}

__host__ void update_values(int* switch_param_host, real* diffu_arr){
	*switch_param_host = 1 - *switch_param_host;
	cudaMemset(diffu_arr, 0, sizeof(real)*MAX_CELL_NUM);
}

void calc_ca2p(CellManager*cm,const real*const ext_stim,const int*const cmap1,const float*const cmap2,real zmax)
{
	//using namespace cont;
	const int offset = cm->nmemb_host + cm->nder_host;
	init_ca2p_avg<<<(cm->ncell_host-offset-1)/128+1,128>>>(cm->ncell_host, offset, cm->state, cm->ca2p_avg);


	const int iz_bound = (int)((zmax + 2.0f*R_max) *inv_dz)+1;
	const int shmem_border = iz_bound + 1 + (iz_bound + 1) % 2;
	static device_alloc_ctor<real> air_stim_arr((NX + 1)*(NY + 1)*(NZ + 1));
	static device_alloc_ctor<real> ATP1((NX + 1)*(NY + 1)*(NZ + 1));
	static device_alloc_ctor<real> ATP2((NX + 1)*(NY + 1)*(NZ + 1));
	static device_alloc_ctor<const real*> cell_diffu_map((NX + 1)*(NY + 1)*(NZ + 1));
	static device_alloc_ctor<real> cell_diffu(MAX_CELL_NUM);
	static device_alloc_ctor<real> dummy_diffu(1);
	static int current_ca2p_switch = 0;
	static real* ATP[2] = { ATP1.ptr, ATP2.ptr };
	cudaMemset(dummy_diffu.ptr, 0, sizeof(real));

	//static RawArr3D<uint_fast8_t> air_stim_flg;
	//static RawArr3D<const double*> cell_diffu_map;
	//const double dummy_diffu = 0;
	init_ca2p_map<<<dim3(NX,NY),iz_bound>>>(air_stim_arr.ptr, cell_diffu_map.ptr, cell_diffu.ptr, cmap1, cm->connection_data, cm->state, iz_bound, dummy_diffu.ptr);
	printf("1 error code:%d\n", cudaGetLastError());
	//test
	cudaError_t ee;
	for (size_t cstp = 0; cstp < Ca_ITR; cstp++) {

		dead_IP3_calc << <(cm->ncell_host - offset - 1) / 128 + 1, 128 >> >(
			cm->ncell_host,
			offset,
			cm->connection_data,
			cm->state,
			cm->IP3[current_ca2p_switch],
			cm->IP3[1 - current_ca2p_switch]);

		
		supra_calc << <(cm->ncell_host - offset - 1) / 128 + 1, 128 >> >(cm->ncell_host, offset, cm->connection_data, cm->current_pos_host(), cm->state, cm->agek, cm->IP3[current_ca2p_switch],
			cm->IP3[1 - current_ca2p_switch], cm->ca2p[current_ca2p_switch], cm->ca2p[1 - current_ca2p_switch], cm->ca2p_avg, cm->ex_inert, cell_diffu.ptr, ext_stim, ATP1.ptr);


		cudaDeviceSynchronize();
		ATP_refresh << <dim3(NX, NY), iz_bound, (shmem_border)*sizeof(int) + (iz_bound + 1)*sizeof(float) >> >(ATP[current_ca2p_switch], ATP[1 - current_ca2p_switch], cell_diffu_map.ptr, air_stim_arr.ptr, cmap2, iz_bound, shmem_border);
		cudaDeviceSynchronize();

		//update_values(cman, ATP);
		update_values(&current_ca2p_switch, cell_diffu.ptr);


	}
	set_cell_ca2p << <(cm->ncell_host - offset - 1) / 128 + 1, 128 >> >(cm->ncell_host, offset, cm->ca2p_avg, cm->state, cm->ca2p[current_ca2p_switch], cm->ca2p[1 - current_ca2p_switch]);
	cudaDeviceSynchronize();
	//system("pause");
	//set_cell_ca2p(cman);
}