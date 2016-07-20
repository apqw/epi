/*
 * CellManager.h
 *
 *  Created on: 2016/07/07
 *      Author: yasu7890v
 */

#ifndef CELLMANAGER_H_
#define CELLMANAGER_H_
#include "define.h"
#include "Lockfree_dev_queue.h"
#include "ROHashmap.h"
#define NEED_RECONNECT (1)
#define NO_NEED_RECONNECT (0)
#include <curand_kernel.h>
class CellDeviceWrapper;

struct CellConnectionData{
	unsigned int connect_num;
	unsigned int last_connect_num;
	ROHashmap<int, MAX_CONNECT_CELL_NUM> order_mem;
	CellIndex connect_index[MAX_CONNECT_CELL_NUM];
	real gj_alloc[2*MAX_CONNECT_CELL_NUM];
	int gj_switch;
	__device__ void add_atomic(CellIndex idx);
	__device__ real* current_gj(){
		return gj_alloc + (gj_switch)*MAX_CONNECT_CELL_NUM;
	}
	__device__ real* last_gj(){
		return gj_alloc + (1-gj_switch)*MAX_CONNECT_CELL_NUM;
	}
	__device__ void switch_gj(){
		gj_switch = 1 - gj_switch;
	}
	CellConnectionData();
};

struct CellDataStruct{

	int* __block_max_store;
	int* current_phase;
	int* need_reconnect;
	//__global__ void set_connection_data_ref

	CELL_STATE* state;
	CellPos* pos[2];
	CellIndex* fix_origin;
	CellConnectionData* connection_data;
	CellIndex* pair_index;
	RealArr ca2p[2];
	RealArr ca2p_avg;
	RealArr ex_inert;
	RealArr IP3[2];
	RealArr agek;
	RealArr ageb;
	RealArr ex_fat;
	RealArr in_fat;
	RealArr spr_nat_len;
	int* uid;
	int* rest_div_times;
	CellIndex* dermis_index;
	real* zmax;
	int* ncell;
	int* nder;
	int* nmemb;

	int* sw;
	int* next_uid;

	static const int DEDUCT_ZMAX_DIV_NUM = 8;
};


class CellManager_Device:public CellDataStruct{
	friend class CellManager;
	void __alloc();
public:
	LockfreeDeviceQueue<CellIndex, MAX_CELL_NUM>* remove_queue;
	__device__ CellPos* current_pos(){
		return pos[*current_phase];
	}
	__device__ CellDeviceWrapper alloc_new_cell();
	__device__ void migrate(CellIndex src, CellIndex dest);
	__device__ void set_need_reconnect(int need);
	
};

class CellManager:public CellDataStruct{

	void alloc();
	void dealloc();
	void memb_init(int,CellConnectionData[]);
	unsigned int* dev_remove_queue_size;
public:
	int current_phase_host;
	int ncell_host;
	int nder_host;
	int nmemb_host;
	int sw_host;
	
	__host__ void set_cell_nums(int _ncell,int _nmemb,int _nder);
		__host__ void fetch_cell_nums();
		__host__ void switch_phase();
	void init_with_file(const char* filename,bool resume=false);
	CellPos* current_pos_host();
	CellPos* next_pos_host();
	CellManager_Device* dev_ptr;
	curandState* rng_state;
	bool should_force_reconnect();
	void no_need_reconnect();
	__host__ unsigned int get_device_remove_queue_size();
	__host__ void reset_device_remove_queue();
	CellManager():current_phase_host(0),ncell_host(0),nder_host(0),nmemb_host(0){
		alloc();
	}
	~CellManager(){
		dealloc();
	}
	
};

//dont keep using over phase
class CellDeviceWrapper{
	CellManager_Device* ptr;
	
	//int current_phase_cache;
	real ca2p_avg_cache;
public:
	const CellIndex index;
	__device__ real& ageb(){ return ptr->ageb[index]; }
	__device__ real& agek(){ return ptr->agek[index]; }
	__device__ CellPos& pos(){ return ptr->pos[*ptr->current_phase][index]; }
	__device__ CellPos& next_pos(){ return ptr->pos[1 - *ptr->current_phase][index]; }
	__device__ CELL_STATE& state(){ return ptr->state[index]; }//__device__ CELL_STATE& state_ref(){ return ptr->state[index]; }
	__device__ CellIndex& fix_origin(){ return ptr->fix_origin[index]; }
	__device__ real& ca2p(){ return ptr->ca2p[*ptr->current_phase][index]; }
	__device__ real& next_ca2p(){ return ptr->ca2p[1-*ptr->current_phase][index]; }
	__device__ real ca2p_avg()const { return ca2p_avg_cache; }__device__ real& ca2p_avg_ref() { return ptr->ca2p_avg[index]; }
	__device__ real& ex_inert(){ return ptr->ex_inert[index]; }
	__device__ real& IP3(){ return ptr->IP3[*ptr->current_phase][index]; }
	__device__ real& next_IP3(){ return ptr->IP3[1-*ptr->current_phase][index]; }
	__device__ real& ex_fat(){ return ptr->ex_fat[index]; }
	__device__ real& in_fat(){ return ptr->in_fat[index]; }
	__device__ real& spr_nat_len(){ return ptr->spr_nat_len[index]; }
	__device__ int& rest_div_times(){ return ptr->rest_div_times[index]; }
	__device__ CellIndex& dermis_index(){ return ptr->dermis_index[index]; }
	__device__ int& uid(){ return ptr->uid[index]; }
	__device__ CellDeviceWrapper(CellManager_Device* dev_ptr, CellIndex _index) :ptr(dev_ptr), index(_index){
		ca2p_avg_cache = dev_ptr->ca2p_avg[_index];
	}
	__device__ CellIndex& pair_index(){ return ptr->pair_index[index]; }
	__device__ CellConnectionData* connection_data_ptr(){ return &(ptr->connection_data[index]); }

	__device__ void remove(){
		ptr->remove_queue->push_back(index);
	}
};
real get_cell_zmax(CellManager*cm);
#endif /* CELLMANAGER_H_ */
