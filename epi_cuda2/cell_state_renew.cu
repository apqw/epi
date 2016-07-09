
#define NDEBUG

#include "define.h"
#include "CellManager.h"
#include <cstdio>
#include <cassert>
#include <curand_kernel.h>
#include "cell_state_renew.h"
#include "utils.h"
#define stoch_div_time_ratio (0.25f)
#define stoch_corr_coef (1.0f)
#define S2 (0.1f)
#define S0 (0.1f*0.2f) //TODO:���킩��₷������
#define delta_L (0.01f*R_max)

#define printf(x,...)
//#define assert(x)
__device__ inline bool __is_malig(CellDeviceWrapper cell){
	return cell.fix_origin() >= 0 && cell.fix_origin() < MALIGNANT;
}

__device__ inline float genrand_real(){
	int id = blockDim.x*blockIdx.x + threadIdx.x;
	curandState state;
	curand_init((unsigned int)clock64(), id, 0, &state);
	return curand_uniform(&state);
}

/**
*  @param [in] c �v�Z�Ώۂ̍זE�זE
*  @return �����ɂ���ĕ␳���ꂽeps_ks
*/
#define eps_ks (0.10f*0.5f)
#define accel_diff (1.0f)
__device__ inline float weighted_eps_ks(CellDeviceWrapper cell) {

	return __is_malig(cell) ? accel_diff *eps_ks : eps_ks;
	//undef?
}

__device__ inline float weighted_eps_ks_impl(float malig_coef) {
	return interp(1.0f,accel_diff,malig_coef)*eps_ks;
	//return (accel_diff*malig_coef+1.0f*(1.0f-malig_coef)) *eps_ks;
	//undef?
}

#define alpha_k (2.0f)
__device__ inline float agek_ALIVE_const_impl(float ca2p_avg,float malig_coef) {
	//using namespace cont;


	return weighted_eps_ks_impl(malig_coef)*(S0 + alpha_k*min0(ca2p_avg - ca2p_init));
}

__device__ inline float agek_ALIVE_const(CellDeviceWrapper cell) {
	//using namespace cont;

	return weighted_eps_ks(cell)*(S0 + alpha_k*min0(cell.ca2p_avg() - ca2p_init));
}

/**
*  @param [in] c �v�Z�Ώۂ̍זE�זE
*  @return �����ɂ���ĕ␳���ꂽeps_kb
*/
#define accel_div (1.0f)
#define eps_kb (0.12f)
__device__ inline float weighted_eps_kb(CellDeviceWrapper cell) {

	return __is_malig(cell) ? accel_div *eps_kb : eps_kb;
}

__device__ inline float weighted_eps_kb_impl(float malig_coef) {

	return interp(1.0f,accel_div,malig_coef)*eps_kb;
}

//using namespace cont;
#define alpha_b (5.0f)
__device__ inline float ageb_const(CellDeviceWrapper cell) {


	return weighted_eps_kb(cell)*(S2 + alpha_b*min0(cell.ca2p_avg() - ca2p_init));
}

__device__ inline float ageb_const_impl(float ca2p_avg,float malig_coef) {

	return weighted_eps_kb_impl(malig_coef)*(S2 + alpha_b*min0(ca2p_avg - ca2p_init));
}


__device__ inline float get_div_age_thresh(CELL_STATE state){
	/** MUSUME�p����J�n�N��̂������l*/
	#define agki_max (6.0f)

	/** FIX�p����J�n�N��̔{�� */
	#define fac (1.0f)

	/** FIX�p����J�n�N��̂������l */
	#define agki_max_fix (fac*agki_max)

	return state == FIX ? agki_max_fix
		: state == MUSUME ? agki_max
		: 0.0f;
}


/**
*  ���􂷂邩�ǂ������m���I�Ɍ��肷�镔��
*  @return ������s���ꍇtrue�A�����łȂ��Ȃ�false
*  @attention �m���I�ȕ�����s��Ȃ��ꍇ�A���true��Ԃ��B
*/
__device__ inline bool stochastic_div_test(CellDeviceWrapper cell) {
	//using namespace cont;
	if (STOCHASTIC!=1) {
		return true;
	}
	else {
		return genrand_real()*(get_div_age_thresh(cell.state())*stoch_div_time_ratio) <= stoch_corr_coef*DT_Cell*weighted_eps_kb(cell)*S2;
	}
}

__device__ inline bool stochastic_div_test_impl(float _div_age_thresh,float malig_coef) {

		return STOCHASTIC!=1 ||
				(genrand_real()*(_div_age_thresh*stoch_div_time_ratio) <= stoch_corr_coef*DT_Cell*weighted_eps_kb_impl(malig_coef)*S2);

}

/**
����̏������ł��Ă��邩�ǂ���
@attention �m���I�ȕ�����s���ꍇ�A�����ɂ��Ԃ�l�͕ς��B
*/
__device__ inline bool is_divide_ready(CellDeviceWrapper cell) {

	if (cell.pair_index()<0 && (cell.ageb() >= get_div_age_thresh(cell.state())*(1.0 - stoch_div_time_ratio))) {
		return stochastic_div_test(cell);
	}
	else {
		return false;
	}
}

__device__ inline bool is_divide_ready_impl(bool has_pair,float ageb,float _div_age_thresh,float malig_coef) {

	if (!has_pair && (ageb >= _div_age_thresh*(1.0 - stoch_div_time_ratio))) {
		return stochastic_div_test_impl(_div_age_thresh,malig_coef);
	}
	else {
		return false;
	}
}

/**
*  c2����c1�ւ̒P�ʃx�N�g�����v�Z����B
*  @param [in] c1 �זEc1
* 	@param [in] c2 �זEc2
*  @param [out] outx �P�ʃx�N�g����x����
*  @param [out] outy �P�ʃx�N�g����y����
*  @param [out] outz �P�ʃx�N�g����z����
*
*  @attention c1��c2�������זE���w���Ă���ꍇ�A����͖���`�B
*  �܂��Aoutx,outy,outz�̂����ꂩ�������ϐ����w���Ă���ꍇ������͖���`�B
*/
__device__ inline float4 calc_cell_uvec(CellPos c1, CellPos c2) {
	const float nvx = p_diff_x(c1.x, c2.x);
	const float nvy = p_diff_y(c1.y, c2.y);
	const float nvz = c1.z - c2.z;

	const float invnorm = rsqrtf(nvx*nvx + nvy*nvy + nvz*nvz);
	return make_float4(
		nvx *invnorm,
		nvy *invnorm,
		nvz *invnorm,0.0f
		);

}

/**
*  �ł��߂�dermis�ɑ΂��Ă̕������(�P�ʃx�N�g��)�������_���Ɍv�Z����B
*
*  @param [in] me ������s���זE
*  @param [in] dermis �ŋߖT��dermis
*  @param [out] outx �P�ʃx�N�g����x����
*  @param [out] outy �P�ʃx�N�g����y����
*  @param [out] outz �P�ʃx�N�g����z����
*
*  @attention me��dermis�������זE���w���Ă���ꍇ�A����͖���`�B
*  �܂��Aoutx,outy,outz�̂����ꂩ�������ϐ����w���Ă���ꍇ������͖���`�B
*/
__device__ float4 div_direction(CellPos c1, CellPos dermis1) {

	//double nx, ny, nz;
	const float4 dermis_normal = calc_cell_uvec(c1, dermis1);
	float ox, oy, oz;
	float sum;
	do {
		//const float rand_theta = M_PI_F*genrand_real();
		float sr1, cr1;
		sincosf(M_PI_F*genrand_real(), &sr1, &cr1);
		float sr2, cr2;
		sincosf(2.0f * M_PI_F*genrand_real(), &sr2, &cr2);
		const float rvx = sr1*cr2;
		const float rvy = sr1*sr2;
		const float rvz = cr1;

		ox = dermis_normal.y*rvz - dermis_normal.z*rvy;
		oy = dermis_normal.z*rvx - dermis_normal.x*rvz;
		oz = dermis_normal.x*rvy - dermis_normal.y*rvx;
	} while ((sum = ox*ox + oy*oy + oz*oz) < 1.0e-14);
	const float invsum = rsqrtf(sum);
	return make_float4(
		ox *invsum,
		oy *invsum,
		oz *invsum,0.0f
		);
}

__device__ inline void cell_divide(CellManager_Device* cmd, CellDeviceWrapper cell) {
	CellDeviceWrapper new_cell = cmd->alloc_new_cell();

	new_cell.state() = cell.state();
	new_cell.fix_origin() = cell.fix_origin();
	new_cell.pos() = cell.pos();
	new_cell.next_pos() = cell.next_pos();
	new_cell.ca2p() = cell.ca2p();
	new_cell.next_ca2p() = cell.ca2p();
	new_cell.ca2p_avg_ref() = cell.ca2p();
	new_cell.IP3() = cell.IP3();
	new_cell.ex_inert() = cell.ex_inert();
	new_cell.agek() = 0.0f;
	new_cell.ageb() = 0.0f; cell.ageb() = 0.0f;
	new_cell.ex_fat() = 0.0f;
	new_cell.in_fat() = 0.0f;
	new_cell.spr_nat_len() = delta_L; cell.spr_nat_len() = delta_L;
	new_cell.rest_div_times() = cell.rest_div_times();
	if (cell.state() == MUSUME){
		new_cell.rest_div_times()--;
		cell.rest_div_times()--;
	}
	new_cell.pair_index() = cell.index;
	cell.pair_index() = new_cell.index;
	if (cell.dermis_index() < 0){
		printf("dermis not found [in divide phase]\n");
		assert(false);
	}
	const float4 div_dir = div_direction(cell.pos(), cmd->current_pos()[cell.dermis_index()]);
	const CellPos my_pos = cell.pos() + (0.5f*delta_L)*div_dir;
	const CellPos pair_pos = new_cell.pos() - (0.5f*delta_L)*div_dir;
	cell.pos() = my_pos; cell.next_pos() = my_pos;
	new_cell.pos() = pair_pos; new_cell.next_pos() = pair_pos;
	printf("new cell detected. cell_index:%d\n", new_cell.index);
}

/**
*  ���זE�̏�ԍX�V�B
*  �����A����A�זE���̍X�V���s���B
*/
__device__ void _MUSUME_state_renew(CellManager_Device* cmd, CellDeviceWrapper cell) {
	if (cell.dermis_index()<0 && cell.pair_index()<0) {
		if (SYSTEM == WHOLE) {
			//printf("ALIVE detected\n");
			cell.state() = ALIVE;
			//musume->cst.k_aging_start_timestep = cman.current_timestep;
		}
		else if (SYSTEM == BASAL) {
			//musume->state = DISA;
			cell.remove();
		}
		return;
	}

	if (cell.rest_div_times() > 0 && is_divide_ready(cell)) {
		cell_divide(cmd, cell);
	}
	else {
		cell.ageb() += DT_Cell*ageb_const(cell);
	}
}
__device__ void _FIX_state_renew(CellManager_Device* cmd, CellDeviceWrapper cell){
if (cell.dermis_index()<0) {
	printf("err: not found fix dermis\n");
	assert(false);
	return;
}
if (is_divide_ready(cell)) {
	cell_divide(cmd, cell);
}
else {
	cell.ageb() += DT_Cell*ageb_const(cell);
}
}

__device__ inline float agek_DEAD_AIR_const() {
	#define eps_kk (0.10f*0.5f)
	#define S1 (0.1f)

	return eps_kk*S1;
}

/**
*  �p�w�A��C�̏�ԍX�V�B
*  ����A�������s���B
*/
__device__ void _DEAD_AIR_state_renew(CellDeviceWrapper cell) {
	#define ADHE_CONST (31.3f)
	#define DISA_conn_num_thresh (11) //Nc

	if (cell.agek() >= ADHE_CONST&&cell.connection_data_ptr()->connect_num <= DISA_conn_num_thresh) {
		cell.remove();
	}
	else {
		cell.agek() += DT_Cell*agek_DEAD_AIR_const();
	}
}

__device__ void cornificate(CellManager_Device* cmd, CellDeviceWrapper cell)
{
	cell.state() = DEAD;
	printf("sw updated:%d\n", atomicAdd(cmd->sw, 1));
}

/** �����̐����E��o��؂�ւ���J���V�E���Z�x */
#define ubar (0.45f)

/**
*  agek�ɂ��X�C�b�`���O�̊ɂ�
*  @note tanh�̕���
*/
#define delta_lipid (0.1f)

/**
*  �J���V�E���Z�x�ɂ��X�C�b�`���O�̊ɂ�
*  @note tanh�̕���
*/
#define delta_sig_r1 (0.1f)

/**
*  @param [in] c �v�Z�Ώۂ̍זE�זE
*  @return ��o���鎉���̎��ԍ���
*/
__device__ inline float k_lipid_release(float ca2p_avg, float agek) {
	//using namespace cont;

	#define lipid_rel (0.05f*4.0f)
	return 0.25f*lipid_rel*(1.0f + tanhf((ca2p_avg - ubar) / delta_sig_r1))*(1.0f + tanhf((agek - THRESH_SP) / delta_lipid));
}

/**
*  @param [in] c �v�Z�Ώۂ̍זE�זE
*  @return �������鎉���̎��ԍ���
*  @attention ���̒l�����ڐ����ʂɂȂ�킯�ł͂Ȃ��B
*/
__device__ inline float k_lipid(float ca2p_avg, float agek) {

	#define lipid (0.6f)
	return 0.25f*lipid*(1 + tanhf((ubar - ca2p_avg) / delta_sig_r1))*(1.0f + tanhf((agek - THRESH_SP) / delta_lipid));
}

__device__ inline void _ALIVE_state_renew(CellManager_Device* cmd, CellDeviceWrapper cell) {
	//using namespace cont;
	const float my_agek = cell.agek();
	//printf("agekno:%f\n", my_agek);
	if (my_agek >= THRESH_DEAD) {
		cornificate(cmd,cell);
	}
	else {
		const float my_in_fat = cell.in_fat();
		const float tmp = k_lipid_release(cell.ca2p_avg(), my_agek)*my_in_fat;
		cell.in_fat() += DT_Cell*(k_lipid(cell.ca2p_avg(), my_agek)*(1.0f - my_in_fat) - tmp);
		cell.ex_fat() += DT_Cell*tmp;
		cell.agek() += DT_Cell*agek_ALIVE_const(cell); //update last
	}
}

__global__ void remove_exec_impl(CellManager_Device* cmd){
	int del_id = blockDim.x*blockIdx.x + threadIdx.x;
	if (del_id < cmd->remove_queue->size()){
		int src = *cmd->ncell - del_id - 1;
		cmd->migrate(src, cmd->remove_queue->data[del_id]);
	}
}

void remove_exec_finalize_dev(CellManager* cm){
	int ncell;
	cudaMemcpy(&ncell, cm->ncell, sizeof(int), cudaMemcpyDeviceToHost);
	ncell -= cm->get_device_remove_queue_size();
	cudaMemcpy(cm->ncell, &ncell, sizeof(int), cudaMemcpyHostToDevice);
	cm->reset_device_remove_queue();
}

__host__ void remove_exec_finalize_host(CellManager* cm){
	cm->fetch_cell_nums();
}

void remove_exec(CellManager* cm){
	remove_exec_impl << <16, cm->get_device_remove_queue_size() / 16 + 1 >> >(cm->dev_ptr);
	//cm->fetch_cell_nums();
	remove_exec_finalize_dev(cm);
	//remove_exec_finalize_host(cm);

}
/*
void state_renew_finalize(CellManager* cm){
	cm->fetch_cell_nums();
}
*/
__global__ inline void cell_state_renew_impl(int ncell,int offset, CellManager_Device* cmd){
	const CellIndex index = blockDim.x*blockIdx.x + threadIdx.x + offset;
	if (index < ncell){
		CellDeviceWrapper c(cmd, index);
		switch (cmd->state[index]){
		case ALIVE:
			_ALIVE_state_renew(cmd, c);
			break;
		case DEAD:case AIR:
			_DEAD_AIR_state_renew(c);
			break;
		case FIX:
			_FIX_state_renew(cmd, c);
			break;
		case MUSUME:
			_MUSUME_state_renew(cmd, c);
			break;
		default:
			break;
		}
	}
}

__device__ void pair_disperse_impl(CellIndex index,CellManager_Device* cmd) {//cannot restrict due to c->pair->pair (== c)

	#define eps_L (0.14f)
	#define unpair_dist_coef (0.9f)
	//using namespace cont;
	//assert(c->pair != nullptr);
	//assert(c->pair->pair == c);
	const float rad_sum = R_max+R_max;
	const float unpair_th = unpair_dist_coef*rad_sum;
	float distSq = 0;
	const CellIndex pairidx = cmd->pair_index[index];
	if (cmd->spr_nat_len[index] < 2.0f*R_max) {
		cmd->spr_nat_len[index] += DT_Cell*eps_L;

		cmd->spr_nat_len[pairidx] = cmd->spr_nat_len[index];
	}
	else if ((distSq = p_cell_dist_sq(cmd->current_pos()[index], cmd->current_pos()[pairidx]))>unpair_th*unpair_th) {
		cmd->spr_nat_len[index] = 0.0f;
		cmd->spr_nat_len[pairidx] = 0.0f;
		cmd->pair_index[pairidx] = -1;
		cmd->pair_index[index] = -1;
		printf("unpaired. distSq:%lf\n", distSq);
	}
}

__global__ inline void pair_disperse(int ncell, int offset, CellManager_Device* cmd){
	const CellIndex index = blockDim.x*blockIdx.x + threadIdx.x + offset;
	if (index < ncell){
		if (cmd->pair_index[index] >= 0)if (index>cmd->pair_index[index]){
			pair_disperse_impl(index, cmd);
		}
	}
}


__global__ inline void cell_value_renew_impl(int ncell,int offset,CellManager_Device*const cmd){

	/*
	 * ex_fat
	 * in_fat
	 * agek
	 * ageb
	 */
	const CellIndex index = blockDim.x*blockIdx.x + threadIdx.x + offset;
	if(index<ncell){
		float agek=cmd->agek[index];
		float ageb=cmd->ageb[index];
		float ex_fat=cmd->ex_fat[index];
		float in_fat=cmd->in_fat[index];
		const float malig_coef=cmd->fix_origin[index]<MALIGNANT?1.0f:0.0f;
		const float my_ca2p_avg=cmd->ca2p_avg[index];

		switch(cmd->state[index]){
		case ALIVE:
			const float tmp = k_lipid_release(my_ca2p_avg, agek)*in_fat;
			in_fat += DT_Cell*(k_lipid(my_ca2p_avg, agek)*(1.0f - in_fat) - tmp);
			ex_fat += DT_Cell*tmp;
			agek   += DT_Cell*agek_ALIVE_const_impl(my_ca2p_avg,malig_coef);
			break;
		case DEAD:case AIR:
			agek+=DT_Cell*agek_DEAD_AIR_const();
			break;
		case MUSUME:case FIX:
			ageb+=DT_Cell*ageb_const_impl(my_ca2p_avg,malig_coef);
			break;
		default:
			break;
		}
		//__syncthreads();
		cmd->agek[index]=agek;
		cmd->ageb[index]=ageb;
		cmd->ex_fat[index]=ex_fat;
		cmd->in_fat[index]=in_fat;
	}
}
__global__ void cell_live_renew_impl(int ncell,int offset,CellManager_Device*const cmd){
	const CellIndex index = blockDim.x*blockIdx.x + threadIdx.x + offset;
		if(index<ncell){
			const CELL_STATE state=cmd->state[index];
			const float agek=cmd->agek[index];
			const float ageb=cmd->ageb[index];
			const CellIndex pair_index=cmd->pair_index[index];
			const CellIndex dermis_index=cmd->dermis_index[index];
			const int rest_div_times=cmd->rest_div_times[index];
			const int fix_origin = cmd->fix_origin[index];
			const int connect_num=cmd->connection_data[index].connect_num;

			const bool has_pair = pair_index>=0;
			const bool musume_diff = state==MUSUME&&dermis_index<0 && !has_pair;
			//(bool has_pair,float ageb,float _div_age_thresh,float malig_coef)
			const bool divide_ready=(state==MUSUME||state==FIX)&&rest_div_times>0
					&&is_divide_ready_impl(has_pair,ageb,get_div_age_thresh(state),fix_origin<MALIGNANT?1.0f:0.0f);
			const bool d_a_remove = (state==DEAD||state==AIR)&&agek >= ADHE_CONST&&connect_num <= DISA_conn_num_thresh;
			const bool musume_remove = musume_diff&&SYSTEM==BASAL;
			const bool musume_to_alive=musume_diff&&SYSTEM==WHOLE;
			const bool alive_cornif = state==ALIVE&&agek>=THRESH_DEAD;
			CellDeviceWrapper cell(cmd, index);
			if(!musume_diff&&divide_ready){
				cell_divide(cmd, cell);
			}

			if(musume_remove||d_a_remove){
				cell.remove();
			}

			if(musume_to_alive){
				cmd->state[index]=ALIVE;
			}

			if(alive_cornif){
				cornificate(cmd,cell);
			}
		}
}

void cell_state_renew(CellManager* cm){
	const int other_num = cm->nder_host + cm->nmemb_host;
	cell_value_renew_impl << <(cm->ncell_host - other_num) / 256+ 1, 256>> >(cm->ncell_host,other_num , cm->dev_ptr);
	cudaDeviceSynchronize();
	cell_live_renew_impl << <(cm->ncell_host - other_num) / 128+ 1, 128>> >(cm->ncell_host,other_num , cm->dev_ptr);
		cudaDeviceSynchronize();
	remove_exec(cm);

	cm->fetch_cell_nums();
	pair_disperse << <(cm->ncell_host - other_num) / 128 + 1, 128 >> >(cm->ncell_host, other_num, cm->dev_ptr);
	cudaDeviceSynchronize();
}
