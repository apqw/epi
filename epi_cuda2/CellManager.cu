/*
 * CellManager.cu
 *
 *  Created on: 2016/07/07
 *      Author: yasu7890v
 */
#include "CellManager.h"
#include "utils.h"
#include <fstream>
#include <vector>
#include <cassert>
#include <string>
__device__ void CellConnectionData::add_atomic(CellIndex idx){
	connect_index[atomicAdd(&connect_num,1)]=idx;
}
#define DCM(ptr,elem_size) cudaMalloc((void**)&(ptr),elem_size*MAX_CELL_NUM)

void CellManager_Device::__alloc(){
	cudaMalloc((void**)&__block_max_store, sizeof(int)*DEDUCT_ZMAX_DIV_NUM);

	DCM(state, sizeof(CELL_STATE));
	DCM(pos[0], sizeof(CellPos));
	DCM(pos[1], sizeof(CellPos));
	DCM(fix_origin, sizeof(CellIndex));
	DCM(connection_data, sizeof(CellConnectionData));
	DCM(pair_index, sizeof(CellIndex));
	DCM(ca2p[0], sizeof(float));
	DCM(ca2p[1], sizeof(float));
	DCM(ca2p_avg, sizeof(float));
	DCM(ex_inert, sizeof(float));
	DCM(IP3, sizeof(float));
	DCM(agek, sizeof(float));
	DCM(ageb, sizeof(float));
	DCM(ex_fat, sizeof(float));
	DCM(in_fat, sizeof(float));
	DCM(spr_nat_len, sizeof(float));
	DCM(rest_div_times, sizeof(int));
	DCM(dermis_index, sizeof(CellIndex));
	cudaMalloc((void**)&zmax, sizeof(float));
	cudaMalloc((void**)&ncell, sizeof(int));
	cudaMalloc((void**)&nder, sizeof(int));
	cudaMalloc((void**)&nmemb, sizeof(int));
	cudaMalloc((void**)&current_phase, sizeof(int));
	cudaMalloc((void**)&sw, sizeof(int));
	cudaMalloc((void**)&need_reconnect, sizeof(int));
	cudaMalloc((void**)&remove_queue, sizeof(LockfreeDeviceQueue<CellIndex, MAX_CELL_NUM>));

	cudaMemset(zmax, 0, sizeof(float));
	cudaMemset(ncell, 0, sizeof(int));
	cudaMemset(nder, 0, sizeof(int));
	cudaMemset(nmemb, 0, sizeof(int));
	cudaMemset(current_phase, 0, sizeof(int));
	cudaMemset(sw, 0, sizeof(int));
	cudaMemset(need_reconnect, 0, sizeof(int));
	LockfreeDeviceQueue<CellIndex, MAX_CELL_NUM> rmq;
	cudaMemcpy(remove_queue, &rmq, sizeof(LockfreeDeviceQueue<CellIndex, MAX_CELL_NUM>), cudaMemcpyHostToDevice);
}


__device__ void CellManager_Device::set_need_reconnect(int need){
	atomicExch(need_reconnect, need);
}
__device__ CellDeviceWrapper CellManager_Device::alloc_new_cell(){
	set_need_reconnect(NEED_RECONNECT);
	return CellDeviceWrapper(this,atomicAdd(ncell, 1));
}

__device__ void CellManager_Device::migrate(CellIndex src,CellIndex dest){
#define MIG(elem) elem[dest]=elem[src]
	MIG(state);
	MIG(pos[0]);
	MIG(pos[1]);
	MIG(fix_origin);
	MIG(connection_data); //slow? invalidated
	MIG(pair_index);
	MIG(ca2p[0]);
	MIG(ca2p[1]);
	MIG(ca2p_avg);
	MIG(ex_inert);
	MIG(IP3);
	MIG(agek);
	MIG(ageb);
	MIG(ex_fat);
	MIG(in_fat);
	MIG(spr_nat_len);
	MIG(rest_div_times);
	MIG(dermis_index);
#undef MIG
	if (pair_index[dest] >= 0){
		pair_index[pair_index[dest]] = dest;
	}

	set_need_reconnect(NEED_RECONNECT);

}
void CellManager::alloc(){
	CellManager_Device tmp;
	tmp.__alloc();


	state=tmp.state;
	pos[0]=tmp.pos[0];
	pos[1]=tmp.pos[1];
	fix_origin=tmp.fix_origin;
	connection_data=tmp.connection_data;
	pair_index=tmp.pair_index;
	ca2p[0]=tmp.ca2p[0];
	ca2p[1] = tmp.ca2p[1];
	ca2p_avg=tmp.ca2p_avg;
	ex_inert=tmp.ex_inert;
	IP3=tmp.IP3;
	agek=tmp.agek;
	ageb=tmp.ageb;
	ex_fat=tmp.ex_fat;
	in_fat=tmp.in_fat;
	spr_nat_len=tmp.spr_nat_len;
	rest_div_times=tmp.rest_div_times;
	dermis_index=tmp.dermis_index;
	zmax=tmp.zmax;
	ncell=tmp.ncell;
	nder=tmp.nder;
	nmemb=tmp.nmemb;
	current_phase=tmp.current_phase;
	sw=tmp.sw;
	need_reconnect = tmp.need_reconnect;
	dev_remove_queue_size = reinterpret_cast<unsigned int*>(tmp.remove_queue); //head address points to size (due to Standard Layout)

	cudaMalloc((void**)&dev_ptr,sizeof(CellManager_Device));
	cudaMemcpy(dev_ptr,&tmp,sizeof(CellManager_Device),cudaMemcpyHostToDevice);

}

void CellManager::dealloc(){
	//CellManager_Device* ptr=new CellManager_Device();

	//cudaMemcpy(ptr,dev_ptr,sizeof(CellManager_Device),cudaMemcpyDeviceToHost);
	cudaFree(__block_max_store);
	cudaFree(state);
	cudaFree(fix_origin);
	cudaFree(connection_data);
	cudaFree(pair_index);
	cudaFree(pos[0]);
	cudaFree(pos[1]);
	cudaFree(ca2p[0]);
	cudaFree(ca2p[1]);
	cudaFree(ca2p_avg);
	cudaFree(ex_inert);
	cudaFree(IP3);
	cudaFree(agek);
	cudaFree(ageb);
	cudaFree(ex_fat);
	cudaFree(in_fat);
	cudaFree(spr_nat_len);
	cudaFree(rest_div_times);
	cudaFree(dermis_index);
	cudaFree(zmax);
	cudaFree(ncell);
	cudaFree(nder);
	cudaFree(nmemb);
	cudaFree(current_phase);
	cudaFree(sw);

	cudaFree(dev_ptr);
	//delete ptr;

}

__host__ void CellManager::memb_init(int _nmemb,CellConnectionData cs[]){
for(int j=0;j<_nmemb;j++){
	const size_t jj=j%NMX;
	const size_t kk = j/NMX;
	cs[j].connect_index[0]=get_adj_memb_idx<DIR_U>(j);
	cs[j].connect_index[1]=get_adj_memb_idx<DIR_L>(j);
	cs[j].connect_index[2]=get_adj_memb_idx<DIR_B>(j);
	cs[j].connect_index[3]=get_adj_memb_idx<DIR_R>(j);
	assert(cs[j].connect_index[0] >= 0 && cs[j].connect_index[1] >= 0 && cs[j].connect_index[2] >= 0 && cs[j].connect_index[3] >= 0);
}
//2-pass
for (int j = 0; j < _nmemb; j++){

	const int memb_u = cs[j].connect_index[0];
	const int memb_l = cs[j].connect_index[1];
	const int memb_b = cs[j].connect_index[2];
	const int memb_r = cs[j].connect_index[3];

	//auto& memb_lu = cs[j].index[4];
	//auto& memb_bl = cs[j].index[5];
	//auto& memb_rb = cs[j].index[6];
	//auto& memb_ur = cs[j].index[7];

	cs[j].connect_index[4] = cs[memb_l].connect_index[0];
	cs[j].connect_index[5] = cs[memb_b].connect_index[1];
	cs[j].connect_index[6] = cs[memb_r].connect_index[2];
	cs[j].connect_index[7] = cs[memb_u].connect_index[3];

	cs[j].connect_num = 8;
	assert(cs[j].connect_index[4] >= 0 && cs[j].connect_index[5] >= 0 && cs[j].connect_index[6] >= 0 && cs[j].connect_index[7] >= 0);
}

}

__host__ void CellManager::init_with_file(const char* filename){
	std::ifstream dstrm(filename);
	if(!dstrm){
		printf("load failed\n");
		assert(dstrm);
	}
	using namespace std;
	vector<CELL_STATE> hstate(MAX_CELL_NUM);
	vector<CellPos> hpos(MAX_CELL_NUM);
	vector<CellIndex> hfix_origin(MAX_CELL_NUM);
	vector<CellConnectionData> hconnection_data(MAX_CELL_NUM);
	vector<CellIndex> hpair_index(MAX_CELL_NUM);
	vector<float> hca2p(MAX_CELL_NUM);
	vector<float>hca2p_avg(MAX_CELL_NUM);
	vector<float>hex_inert(MAX_CELL_NUM);
	vector<float>hagek(MAX_CELL_NUM);
	vector<float>hageb(MAX_CELL_NUM);
	vector<float>hex_fat(MAX_CELL_NUM);
	vector<float>hin_fat(MAX_CELL_NUM);
	vector<float>hspr_nat_len(MAX_CELL_NUM);
	vector<float>hrad(MAX_CELL_NUM);
	vector<int> hrest_div_times(MAX_CELL_NUM);

	string line;
	int id_count=0;
	unsigned int phase=0;
	int nmemb=0;
	int nder=0;

	while (std::getline(dstrm, line)) {

			sscanf(line.c_str(), "%*d %d %f %f %f %f %f %f %f %f %d %f %f %*d %f %d %d",
				&hstate[id_count],
				&hrad[id_count],
				&hageb[id_count],
				&hagek[id_count],
				&hca2p[id_count],
				&hpos[id_count].x,
				&hpos[id_count].y,
				&hpos[id_count].z,
				&hca2p_avg[id_count],
				&hrest_div_times[id_count],
				&hex_fat[id_count],
				&hin_fat[id_count],
				&hspr_nat_len[id_count],
				&hpair_index[id_count],
				&hfix_origin[id_count]);
			hpos[id_count].w = 0.0f;
			/*
			BLANK�ɂ��ǂ蒅������I��
			*/
			CELL_STATE state = hstate[id_count];
			float rad = hrad[id_count];
			if (state == BLANK)break;

			/*
			validation
			*/
			if (SYSTEM == BASAL && (state == ALIVE || state == DEAD || state == AIR)) {
				printf(" input date must not contain ALIVE or DEAD in case of BASAL\n");
				exit(1);
			}
			if (state == DER && rad != R_der) {
				printf("radii of DER not consistent with param.h %lf %lf\n",rad,R_der);
				exit(1);
			}
			if (state == MEMB && rad != R_memb) {
				printf("radii of DER not consistent with param.h\n");
				exit(1);
			}
			if (phase == 0 && state != MEMB) {
				assert(state == DER);
				phase++;
			}

			if (phase == 1 && state != DER) {
				//assert(state == DER);
				phase++;
			}
			if (phase > 0 && state == MEMB) {
				printf("non phase0 memb\n");
				exit(1);
			}
			if (phase > 1 && state == DER) {
				printf("non phase1 der\n");
				exit(1);
			}

			//if (state == FIX)printf("FIX\n");

			if (state == MEMB)nmemb++;
			if (state == DER)nder++;
			//*(unsigned int*)(&c_pos_h[id_count].w) = c_state_h[id_count];
			printf("Phase %d  Cell loaded:%d\n", phase, id_count++);

		}
	assert(nmemb==NMX*NMY);
	printf("eusyo0\n");
	set_cell_nums(id_count,nmemb,nder);
	printf("eusyo1\n");
	memb_init(nmemb, &hconnection_data[0]);
	cudaMemcpy(state, &hstate[0], sizeof(CELL_STATE)*MAX_CELL_NUM, cudaMemcpyHostToDevice);
		//cudaMemcpy(d->c_radius_d, c_radius_h, sizeof(float)*MAX_CELL_NUM, cudaMemcpyHostToDevice);
		cudaMemcpy(ageb, &hageb[0], sizeof(float)*MAX_CELL_NUM, cudaMemcpyHostToDevice);
		cudaMemcpy(agek, &hagek[0], sizeof(float)*MAX_CELL_NUM, cudaMemcpyHostToDevice);
		cudaMemcpy(ca2p[0], &hca2p[0], sizeof(float)*MAX_CELL_NUM, cudaMemcpyHostToDevice);
		cudaMemcpy(ca2p[1], &hca2p[0], sizeof(float)*MAX_CELL_NUM, cudaMemcpyHostToDevice);
		cudaMemcpy(pos[0], &hpos[0], sizeof(CellPos)*MAX_CELL_NUM, cudaMemcpyHostToDevice);
		cudaMemcpy(pos[1], &hpos[0], sizeof(CellPos)*MAX_CELL_NUM, cudaMemcpyHostToDevice);
		cudaMemcpy(ca2p_avg, &hca2p_avg[0], sizeof(float)*MAX_CELL_NUM, cudaMemcpyHostToDevice);
		cudaMemcpy(rest_div_times, &hrest_div_times[0], sizeof(int)*MAX_CELL_NUM, cudaMemcpyHostToDevice);
		cudaMemcpy(ex_fat, &hex_fat[0], sizeof(float)*MAX_CELL_NUM, cudaMemcpyHostToDevice);
		cudaMemcpy(in_fat, &hin_fat[0], sizeof(float)*MAX_CELL_NUM, cudaMemcpyHostToDevice);
		cudaMemcpy(spr_nat_len, &hspr_nat_len[0], sizeof(float)*MAX_CELL_NUM, cudaMemcpyHostToDevice);
		cudaMemcpy(pair_index, &hpair_index[0], sizeof(int)*MAX_CELL_NUM, cudaMemcpyHostToDevice);
		cudaMemcpy(fix_origin, &hfix_origin[0], sizeof(int)*MAX_CELL_NUM, cudaMemcpyHostToDevice);
		cudaMemcpy(connection_data, &hconnection_data[0], sizeof(CellConnectionData)*MAX_CELL_NUM, cudaMemcpyHostToDevice);
}

__host__ void CellManager::switch_phase(){
	current_phase_host = 1 - current_phase_host;
	cudaMemcpy(current_phase, &current_phase_host, sizeof(int), cudaMemcpyHostToDevice);
}

__host__ void CellManager::fetch_cell_nums(){
	cudaMemcpy(&ncell_host,ncell,sizeof(int),cudaMemcpyDeviceToHost);
	cudaMemcpy(&nder_host,nder,sizeof(int),cudaMemcpyDeviceToHost);
	cudaMemcpy(&nmemb_host,nmemb,sizeof(int),cudaMemcpyDeviceToHost);
}

bool CellManager::should_force_reconnect(){
	int tmp;
	cudaMemcpy(&tmp, need_reconnect, sizeof(int), cudaMemcpyDeviceToHost);
	return ( tmp == NEED_RECONNECT);
}

void CellManager::no_need_reconnect(){
	int tmp = NO_NEED_RECONNECT;
	cudaMemcpy(need_reconnect, &tmp, sizeof(int), cudaMemcpyHostToDevice);
}
unsigned int CellManager::get_device_remove_queue_size(){
	unsigned int tmp;
	cudaMemcpy(&tmp, dev_remove_queue_size, sizeof(unsigned int), cudaMemcpyDeviceToHost);
	return tmp;
}

void CellManager::reset_device_remove_queue(){
	//unsigned int tmp;
	cudaMemset(dev_remove_queue_size,0, sizeof(unsigned int));
	//return tmp;
}
__host__ void CellManager::set_cell_nums(int _ncell,int _nmemb,int _nder){

	ncell_host=_ncell;
	nmemb_host=_nmemb;

	nder_host = _nder;

printf("eruusi2\n");
	cudaMemcpy(ncell,&ncell_host,sizeof(int),cudaMemcpyHostToDevice);
		cudaMemcpy(nder,&nder_host,sizeof(int),cudaMemcpyHostToDevice);
		cudaMemcpy(nmemb,&nmemb_host,sizeof(int),cudaMemcpyHostToDevice);
}
CellPos* CellManager::current_pos_host(){
	return pos[current_phase_host];
}

CellPos* CellManager::next_pos_host(){
	return pos[1-current_phase_host];
}
