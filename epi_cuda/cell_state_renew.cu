#include "global_cuda.h"
#include "utils_cuda.cuh"
#include "remove_queue.cuh"
#include <curand_kernel.h>
#include <cassert>
#include "cell_state_renew.h"
#define stoch_div_time_ratio (0.25f)

//only for 1-d block and 1-d thread
__device__ inline float genrand_real(){
	int id = blockDim.x*blockIdx.x + threadIdx.x;
	curandState state;
	curand_init((unsigned int)clock64(), id, 0, &state);
	return curand_uniform(&state);
}

////////////////////////////////////////////////////////////////////////////////
//                                                                            
//    agekの計算の定義                                                        
//                                                                            
////////////////////////////////////////////////////////////////////////////////



/**
*  @param [in] c 計算対象の細胞細胞
*  @return 条件によって補正されたeps_ks
*/
__device__ inline float weighted_eps_ks(bool is_malignant=false) {
	const float eps_ks = 0.10f*0.5f;
	const float accel_diff = 1.0f;
	return is_malignant ? accel_diff *eps_ks : eps_ks;
}

__device__ inline float agek_ALIVE_const(float ca2p_avg, bool is_malignant=false) {
	using namespace cont;
	const float alpha_k = 2.0f;

	return weighted_eps_ks(is_malignant)*(S0 + alpha_k*min0(ca2p_avg - ca2p_init));
}

__device__ inline float agek_DEAD_AIR_const() {
	const float eps_kk = 0.10f*0.5f;
	const float S1 = 0.1f;

	return eps_kk*S1;
}

////////////////////////////////////////////////////////////////////////////////
//                                                                            
//    agebの計算の定義                                                        
//                                                                            
////////////////////////////////////////////////////////////////////////////////

#define S2 (0.1f) //TODO:よりわかりやすい命名

/**
*  @param [in] c 計算対象の細胞細胞
*  @return 条件によって補正されたeps_kb
*/
__device__ inline float weighted_eps_kb(bool is_malignant=false) {
	const float accel_div = 1.0f;
	const float eps_kb = 0.12f;
	return is_malignant ? accel_div *eps_kb : eps_kb;
}

__device__ inline float ageb_const(float ca2p_avg, bool is_malignant=false) {
	using namespace cont;
	const float alpha_b = 5.0;

	return weighted_eps_kb(is_malignant)*(S2 + alpha_b*min0(ca2p_avg - ca2p_init));
}


////////////////////////////////////////////////////////////////////////////////
//                                                                            
//    細胞分裂に関する定義                                                    
//                                                                            
////////////////////////////////////////////////////////////////////////////////

/**
*  分裂するかどうかを確率的に決定する部分
*  @return 分裂を行う場合true、そうでないならfalse
*  @attention 確率的な分裂を行わない場合、常にtrueを返す。
*/
__device__ inline bool stochastic_div_test(float div_age_thresh,bool is_malignant=false) {
	using namespace cont;
	if (STOCHASTIC != 1) {
		return true;
	}
	else {
		return genrand_real()*(div_age_thresh*stoch_div_time_ratio) <= stoch_corr_coef*DT_Cell*weighted_eps_kb(is_malignant)*S2;
	}
}


__device__ inline float get_div_age_thresh(unsigned int state){
	/** MUSUME用分裂開始年齢のしきい値*/
	const float agki_max = 6.0f;

	/** FIX用分裂開始年齢の倍率 */
	const float fac = 1.0f;

	/** FIX用分裂開始年齢のしきい値 */
	const float agki_max_fix = fac*agki_max;
	return state == FIX ? agki_max_fix
		: state == MUSUME ? agki_max
		: 0.0f;
}
/**
分裂の準備ができているかどうか
@attention 確率的な分裂を行う場合、乱数により返り値は変わる。
*/
__device__ inline bool is_divide_ready(unsigned int state, bool has_pair, float ageb, bool is_malignant = false) {
	const float div_age_thresh = get_div_age_thresh(state);
	if (!has_pair && (ageb >= div_age_thresh*(1.0f - stoch_div_time_ratio))) {
		return stochastic_div_test(div_age_thresh,is_malignant);
	}
	else {
		return false;
	}
}

/**
*  c2からc1への単位ベクトルを計算する。
*  @param [in] c1 細胞c1
* 	@param [in] c2 細胞c2
*  @param [out] outx 単位ベクトルのx成分
*  @param [out] outy 単位ベクトルのy成分
*  @param [out] outz 単位ベクトルのz成分
*
*  @attention c1とc2が同じ細胞を指している場合、動作は未定義。
*  また、outx,outy,outzのいずれかが同じ変数を指している場合も動作は未定義。
*/
__device__ inline float3 calc_cell_uvec(cell_pos_set c1,cell_pos_set c2) {
	const float nvx = p_diff_x(c1.x, c2.x);
	const float nvy = p_diff_y(c1.y, c2.y);
	const float nvz = c1.z - c2.z;

	const float norm = sqrtf(nvx*nvx+nvy*nvy+nvz*nvz);
	return make_float3(
		nvx / norm,
		nvy / norm,
		nvz / norm
		);

}

/**
*  最も近いdermisに対しての分裂方向(単位ベクトル)をランダムに計算する。
*
*  @param [in] me 分裂を行う細胞
*  @param [in] dermis 最近傍のdermis
*  @param [out] outx 単位ベクトルのx成分
*  @param [out] outy 単位ベクトルのy成分
*  @param [out] outz 単位ベクトルのz成分
*
*  @attention meとdermisが同じ細胞を指している場合、動作は未定義。
*  また、outx,outy,outzのいずれかが同じ変数を指している場合も動作は未定義。
*/
__device__ float3 div_direction(cell_pos_set c1,cell_pos_set dermis1) {

	//double nx, ny, nz;
	const float3 dermis_normal=calc_cell_uvec(c1, dermis1);
	float ox, oy, oz;
	float sum;
	do {
		const float rand_theta = M_PI_F*genrand_real();
		const float cr1 = cosf(rand_theta);
		const float sr1 = sinf(rand_theta);
		const float rand_phi = 2.0f * M_PI_F*genrand_real();
		const float cr2 = cosf(rand_phi);
		const float sr2 = sinf(rand_phi);
		const float rvx = sr1*cr2;
		const float rvy = sr1*sr2;
		const float rvz = cr1;

		ox = dermis_normal.y*rvz - dermis_normal.z*rvy;
		oy = dermis_normal.z*rvx - dermis_normal.x*rvz;
		oz = dermis_normal.x*rvy - dermis_normal.y*rvx;
	} while ((sum = ox*ox+oy*oy+oz*oz) < 1.0e-14);
	sum = sqrt(sum);
	return make_float3(
	ox / sum,
	oy / sum,
	oz / sum
	);
}
/**
*  細胞分裂を行う。
*
*  @param cman cellmanager
*  @param div 分裂を行う細胞
*
*  @attention 最近傍のdermisが計算されている必要がある。
*/
__device__ inline void cell_divide(
	int index,
	int current_switch,
	int dermis,
	int* ncell,
	unsigned int cstate[],
	int fix_origin[],
	cell_pos_set* pos_set[],
	float* ca2p[],
	float ca2p_avg[],
	float ex_inert[],
	float IP3[],
	float agek[],
	float ageb[],
	float ex_fat[],
	float in_fat[],
	float spr_nat_len[],
	int rest_div_times[],
	int pair_index[]) {


	using namespace cont;
	const float delta_L = 0.01f*R_max;
	int nc_idx = atomicAdd(ncell, 1);
	cstate[nc_idx] = MUSUME;
	fix_origin[nc_idx] = fix_origin[index];
	pos_set[current_switch][nc_idx] = pos_set[current_switch][index];
	pos_set[1 - current_switch][nc_idx] = pos_set[1 - current_switch][index];
	//no rad
	ca2p[0][nc_idx] = ca2p[current_switch][index];
	ca2p[1][nc_idx] = ca2p[current_switch][index]; //0
	ca2p_avg[nc_idx] = ca2p[current_switch][index];
	IP3[nc_idx] = IP3[index];
	ex_inert[nc_idx] = ex_inert[index];
	agek[nc_idx] = 0.0f;
	ageb[nc_idx] = 0.0f; ageb[index] = 0.0f;
	ex_fat[nc_idx] = 0.0f;
	in_fat[nc_idx] = 0.0f;
	spr_nat_len[nc_idx] = delta_L; spr_nat_len[index] = 0.0f;
	rest_div_times[nc_idx] = rest_div_times[index];
	if (cstate[index] == MUSUME){
		rest_div_times[nc_idx]--;
		rest_div_times[index]--;
	}
	pair_index[nc_idx] = index;
	pair_index[index] = nc_idx;
	if (dermis < 0){
		printf("dermis not found [in divide phase]\n");
		assert(dermis >= 0);
	}
	/*
	div->pair = cman.create(
		MUSUME,
		div->fix_origin,
		div->x(), div->y(), div->z(),
		div->radius,
		div->ca2p(),
		div->ca2p(),//this is avg value,do not use orignal avg
		div->IP3(),
		div->ex_inert,
		0, 0,//set ages 0
		0, 0,//set fats 0
		delta_L,//init
		div->rest_div_times,
		div->is_malignant
		);
		*/
	/*
	div->pair->pair = div;
	div->spr_nat_len = delta_L;
	div->ageb = 0;
	if (div->state == MUSUME) {
		div->rest_div_times--;
		div->pair->rest_div_times--;
	}
	if (div->dermis() == nullptr) {
		std::cout << "No dermis found in divide_try." << std::endl;
		exit(1);
	}
	*/
	//assert(div->dermis() != nullptr);
	float divx, divy, divz;
	const float3 div_dir = div_direction(pos_set[current_switch][index], pos_set[current_switch][dermis]);
	//!set value
	const cell_pos_set my_pos = make_float4(
		pos_set[current_switch][index].x + div_dir.x*0.5f*delta_L,
		pos_set[current_switch][index].y + div_dir.y*0.5f*delta_L,
		pos_set[current_switch][index].z + div_dir.z*0.5f*delta_L,
		pos_set[current_switch][index].w //state
		);

	const cell_pos_set pair_pos = make_float4(
		pos_set[current_switch][nc_idx].x - div_dir.x*0.5f*delta_L,
		pos_set[current_switch][nc_idx].y - div_dir.y*0.5f*delta_L,
		pos_set[current_switch][nc_idx].z - div_dir.z*0.5f*delta_L,
		0.0f
		);
	*(unsigned int*)&pair_pos.w = MUSUME;


	pos_set[current_switch][index] = my_pos;
	pos_set[1 - current_switch][index] = my_pos;
	pos_set[current_switch][nc_idx] = pair_pos;
	pos_set[1 - current_switch][nc_idx] = pair_pos;
	/*
	+=
	div->x._set(div->x() + divx*0.5*delta_L);
	div->y._set(div->y() + divy*0.5*delta_L);
	div->z._set(div->z() + divz*0.5*delta_L);

	div->pair->x._set(div->pair->x() - divx*0.5*delta_L);
	div->pair->y._set(div->pair->y() - divy*0.5*delta_L);
	div->pair->z._set(div->pair->z() - divz*0.5*delta_L);
	*/
	printf("new cell detected. cell_num:%zd\n", nc_idx+1);
}

////////////////////////////////////////////////////////////////////////////////
//                                                                            
//    脂質の計算の定義                                                        
//                                                                            
////////////////////////////////////////////////////////////////////////////////

/** 脂質の生成・放出を切り替えるカルシウム濃度 */
#define ubar (0.45f)

/**
*  agekによるスイッチングの緩さ
*  @note tanhの分母
*/
#define delta_lipid (0.1f)

/**
*  カルシウム濃度によるスイッチングの緩さ
*  @note tanhの分母
*/
#define delta_sig_r1 (0.1f)

/**
*  @param [in] c 計算対象の細胞細胞
*  @return 放出する脂質の時間差分
*/
__device__ inline float k_lipid_release(float ca2p_avg, float agek) {
	using namespace cont;

	const float lipid_rel = 0.05f*4.0f;
	return 0.25f*lipid_rel*(1.0f + tanhf((ca2p_avg - ubar) / delta_sig_r1))*(1.0f + tanh((agek - THRESH_SP) / delta_lipid));
}

/**
*  @param [in] c 計算対象の細胞細胞
*  @return 生成する脂質の時間差分
*  @attention この値が直接生成量になるわけではない。
*/
__device__ inline float k_lipid(float ca2p_avg,float agek) {
	using namespace cont;

	const float lipid = 0.6f;
	return 0.25f*lipid*(1 + tanhf((ubar - ca2p_avg) / delta_sig_r1))*(1.0f + tanh((agek - THRESH_SP) / delta_lipid));
}

__device__ void cornificate(int index,unsigned int* cstate,int* sw)
{
	cstate[index] = DEAD;
	//al->cst.k_cornified_timestep = cman.current_timestep;
	printf("sw updated:%d\n", atomicAdd(sw,1));
}

__device__ inline void _ALIVE_state_renew(int index,bool is_malignant,float ca2p_avg,unsigned int* cstate,float* agek,float* in_fat,float* ex_fat,int* sw) {
	using namespace cont;
	const float my_agek = agek[index];
	//printf("agekno:%f\n", my_agek);
	if (my_agek >= THRESH_DEAD) {
		cornificate(index,cstate,sw);
	}
	else {
		const float my_in_fat = in_fat[index];
		const float tmp = k_lipid_release(ca2p_avg, my_agek)*my_in_fat;
		in_fat[index] += DT_Cell*(k_lipid(ca2p_avg, my_agek)*(1.0f - my_in_fat) - tmp);
		ex_fat[index] += DT_Cell*tmp;
		agek[index] += DT_Cell*agek_ALIVE_const(is_malignant,ca2p_avg); //update last
	}
}

/**
*  角層、空気の状態更新。
*  加齢、剥離を行う。
*/
__device__ void _DEAD_AIR_state_renew(int index, float* agek,unsigned int connected_cell_num,remove_queue* rmq) {
	const float ADHE_CONST = 31.3f;
	const unsigned int DISA_conn_num_thresh = 11; //Nc

	if (agek[index] >= ADHE_CONST&&connected_cell_num <= DISA_conn_num_thresh) {
		//da->state = DISA;
		//printf("eusyo\n");
		rmq->push_back(index);
	}
	else {
		agek[index] += cont::DT_Cell*agek_DEAD_AIR_const();
	}
}

/**
*  幹細胞の状態更新。
*  分裂、細胞周期の更新を行う。
なんとかする
*/
__device__ void _FIX_state_renew(int index,
	int current_switch,
	int dermis_idx,
	int* ncell,
	unsigned int cstate[],
	int fix_origin[],
	cell_pos_set* pos_set[],
	float* ca2p[],
	float ca2p_avg[],
	float ex_inert[],
	float IP3[],
	float agek[],
	float ageb[],
	float ex_fat[],
	float in_fat[],
	float spr_nat_len[],
	int rest_div_times[],
	int pair_index[]) {
	if (dermis_idx<0) {
		printf("err: not found fix dermis\n");
		//printf("x:%lf,y:%lf,z:%lf\n", fix->x(), fix->y(), fix->z());
		//printf("connected_num:%zd\n", fix->connected_cell.size());
		assert(dermis_idx>=0);
		//exit(1);
		return;
	}
	bool is_malignant = fix_origin[index] < MALIGNANT;
	if (is_divide_ready(FIX,pair_index[index]>=0,ageb[index],is_malignant)) {
		cell_divide(index,
			current_switch,
			dermis_idx,
			ncell,
			 cstate,
			fix_origin,
			pos_set,
			ca2p,
			ca2p_avg,
			ex_inert,
			IP3,
			agek,
			ageb,
			ex_fat,
			in_fat,
			spr_nat_len,
			rest_div_times,
			pair_index);
	}
	else {
		ageb[index] += cont::DT_Cell*ageb_const(ca2p_avg[index], is_malignant);
	}
}

/**
*  娘細胞の状態更新。
*  分化、分裂、細胞周期の更新を行う。
*/
__device__  void _MUSUME_state_renew(int index,
	remove_queue* rmq,
	int current_switch,
	int dermis_idx,
	int* ncell,
	unsigned int cstate[],
	int fix_origin[],
	cell_pos_set* pos_set[],
	float* ca2p[],
	float ca2p_avg[],
	float ex_inert[],
	float IP3[],
	float agek[],
	float ageb[],
	float ex_fat[],
	float in_fat[],
	float spr_nat_len[],
	int rest_div_times[],
	int pair_index[]) {
	if (dermis_idx<0 && pair_index[index]<0) {
		if (SYSTEM == WHOLE) {
			printf("ALIVE detected\n");
			cstate[index] = ALIVE;
			*(unsigned int*)&pos_set[0][index].w = ALIVE;
			*(unsigned int*)&pos_set[1][index].w = ALIVE;
			//musume->cst.k_aging_start_timestep = cman.current_timestep;
		}
		else if (SYSTEM == BASAL) {
			//musume->state = DISA;
			rmq->push_back(index);
		}
		return;
	}
	bool is_malignant = fix_origin[index] < MALIGNANT;
	if (rest_div_times[index] > 0 && is_divide_ready(MUSUME, pair_index[index] >= 0, ageb[index], is_malignant)) {
		cell_divide(index,
			current_switch,
			dermis_idx,
			ncell,
			cstate,
			fix_origin,
			pos_set,
			ca2p,
			ca2p_avg,
			ex_inert,
			IP3,
			agek,
			ageb,
			ex_fat,
			in_fat,
			spr_nat_len,
			rest_div_times,
			pair_index);
	}
	else {
		ageb[index] += cont::DT_Cell*ageb_const(ca2p_avg[index], is_malignant);
	}
}
__global__ void remove_exec(remove_queue* rmq,
	int current_ncell,
	int* ncell,
	unsigned int* c_state_d,
	cell_pos_set* c_pos_d[2], //w:=state(duplicated)
	int* c_fix_origin_d,
	connected_index_set* c_connected_index_d,
	int* c_pair_index_d,
	float* c_ca2p_d[2],
	float* c_ca2p_avg_d,
	float* c_ex_inert_d,
	float* c_IP3_d,
	float* c_agek_d,
	float* c_ageb_d,
	float* c_ex_fat_d,
	float* c_in_fat_d,
	float* c_spr_nat_len_d,
	//float radius[];
	int* c_rest_div_times_d,
	int* c_dermis_index_d){
	for (int i = 0; i < rmq->queue_head(); i++){
		int idx = rmq->queue[i];
		int last = current_ncell - i - 1;
#define MIGRATE(arr,rm,last) arr[rm]=arr[last]
		
		MIGRATE(c_state_d, idx, last);
		/*
		MIGRATE(c_pos_d[0], idx, last);
		MIGRATE(c_pos_d[1], idx, last);
		MIGRATE(c_fix_origin_d, idx, last);
		MIGRATE(c_connected_index_d, idx, last); //heavy??
		MIGRATE(c_pair_index_d, idx, last);
		MIGRATE(c_ca2p_d[0], idx, last);
		MIGRATE(c_ca2p_d[1], idx, last);
		MIGRATE(c_ca2p_avg_d, idx, last);
		MIGRATE(c_ex_inert_d, idx, last);
		MIGRATE(c_IP3_d, idx, last);
		MIGRATE(c_agek_d, idx, last);
		MIGRATE(c_ageb_d, idx, last);
		MIGRATE(c_ex_fat_d, idx, last);
		MIGRATE(c_in_fat_d, idx, last);
		MIGRATE(c_spr_nat_len_d, idx, last);
		MIGRATE(c_rest_div_times_d, idx, last);
		MIGRATE(c_dermis_index_d, idx, last);
		*/
	}
	*ncell -= rmq->queue_head();
	rmq->reset();
}
__global__ void cell_state_renew_impl(
	int offset,
	remove_queue* rmq,
	int current_switch,
	int* dermis_idx,
	int current_ncell,
	int* ncell,
	unsigned int cstate[],
	int fix_origin[],
	cell_pos_set* pos_set[],
	float* ca2p[],
	float ca2p_avg[],
	float ex_inert[],
	float IP3[],
	float agek[],
	float ageb[],
	float ex_fat[],
	float in_fat[],
	float spr_nat_len[],
	int rest_div_times[],
	int pair_index[],
	connected_index_set* conn,
	int* sw)
{
	int index = blockDim.x*blockIdx.x + threadIdx.x+offset;
	
	if (index < current_ncell){
		unsigned int state = cstate[index];
		switch (state){
		case ALIVE:
			_ALIVE_state_renew(index, fix_origin[index] < MALIGNANT, ca2p_avg[index], cstate, agek, in_fat, ex_fat, sw);
			break;
		case DEAD:case AIR:
			_DEAD_AIR_state_renew(index, agek, conn[index].connected_num, rmq);
			break;
		case FIX:
			
			_FIX_state_renew(index,
				current_switch,
				dermis_idx[index],
				ncell,
				cstate,
				fix_origin,
				pos_set,
				ca2p,
				ca2p_avg,
				ex_inert,
				IP3,
				agek,
				ageb,
				ex_fat,
				in_fat,
				spr_nat_len,
				rest_div_times,
				pair_index);
				
			break;
		case MUSUME:
			/*
			_MUSUME_state_renew(index,
				rmq,
				current_switch,
				dermis_idx[index],
				ncell,
				cstate,
				fix_origin,
				pos_set,
				ca2p,
				ca2p_avg,
				ex_inert,
				IP3,
				agek,
				ageb,
				ex_fat,
				in_fat,
				spr_nat_len,
				rest_div_times,
				pair_index);
				*/
			break;
		default:
			break;
		}
	}
	/*
	remove_exec(rmq,
		ncell,
		cstate,
		pos_set, //w:=state(duplicated)
		fix_origin,
		conn,
		pair_index,
		ca2p,
		ca2p_avg,
		ex_inert,
		IP3,
		agek,
		ageb,
		ex_fat,
		in_fat,
		spr_nat_len,
		//float radius[];
		rest_div_times,
		dermis_idx);
		*/
	//cman.remove_exec();
	/*
	cman.other_foreach_parallel_native([](Cell*const RESTRICT c) {

		if (c->pair != nullptr) if (no_double_count(c, c->pair)) {
			pair_disperse(c);
		}

	});
	*/
}

void cell_state_renew(DeviceData*d,remove_queue* rmq){
	/*
		remove_queue* rmq,
	int current_switch,
	int* dermis_idx,
	int current_ncell,
	int* ncell,
	unsigned int cstate[],
	int fix_origin[],
	cell_pos_set* pos_set[],
	float* ca2p[],
	float ca2p_avg[],
	float ex_inert[],
	float IP3[],
	float agek[],
	float ageb[],
	float ex_fat[],
	float in_fat[],
	float spr_nat_len[],
	int rest_div_times[],
	int pair_index[],
	connected_index_set* conn,
	int* sw
	*/
	// << <(d->ncell - d->nmemb - d->nder) / 256 + 1, 256 >> >
	int* changed_ncell;
	cudaMalloc(&changed_ncell, sizeof(int));
	//int tmp = 0;
	cudaMemset(changed_ncell, 0, sizeof(int));
	cell_state_renew_impl << <(d->ncell - d->nmemb - d->nder) / 256 + 1, 256 >> >
		(d->nmemb+d->nder,rmq,
		d->current,
		d->c_dermis_index_d,
		d->ncell,
		changed_ncell,
		(unsigned int*)d->c_state_d,
		d->c_fix_origin_d,
		d->c_pos_d,
		d->c_ca2p_d,
		d->c_ca2p_avg_d,
		d->c_ex_inert_d,
		d->c_IP3_d,
		d->c_agek_d,
		d->c_ageb_d,
		d->c_ex_fat_d,
		d->c_in_fat_d,
		d->c_spr_nat_len_d,
		d->c_rest_div_times_d,
		d->c_pair_index_d,
		d->c_connected_index_d,
		d->sw
		);
	cudaThreadSynchronize();
	//has bug
	
	remove_exec << <1, 1 >> >(
		rmq,
		d->ncell,
		changed_ncell,
		(unsigned int*)d->c_state_d,
		d->c_pos_d,
		d->c_fix_origin_d,
		d->c_connected_index_d,
		d->c_pair_index_d,
		d->c_ca2p_d,
		d->c_ca2p_avg_d,
		d->c_ex_inert_d,
		d->c_IP3_d,
		d->c_agek_d,
		d->c_ageb_d,
		d->c_ex_fat_d,
		d->c_in_fat_d,
		d->c_spr_nat_len_d,
		d->c_rest_div_times_d,
		d->c_dermis_index_d
		);
		
	int tmp;
	cudaMemcpy(&tmp, changed_ncell, sizeof(int), cudaMemcpyDeviceToHost);
	d->ncell += tmp;

}