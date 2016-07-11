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
#define NEED_RECONNECT (1)
#define NO_NEED_RECONNECT (0)
#include <curand_kernel.h>
class CellDeviceWrapper;
struct CellConnectionData{
	unsigned int connect_num;
	CellIndex connect_index[MAX_CONNECT_CELL_NUM];
	float gj[MAX_CONNECT_CELL_NUM];
	__device__ void add_atomic(CellIndex idx);
	CellConnectionData():connect_num(0){}
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
	FloatArr ca2p[2];
	FloatArr ca2p_avg;
	FloatArr ex_inert;
	FloatArr IP3[2];
	FloatArr agek;
	FloatArr ageb;
	FloatArr ex_fat;
	FloatArr in_fat;
	FloatArr spr_nat_len;
	int* rest_div_times;
	CellIndex* dermis_index;
	float* zmax;
	int* ncell;
	int* nder;
	int* nmemb;

	int* sw;

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
	
	__host__ void set_cell_nums(int _ncell,int _nmemb,int _nder);
		__host__ void fetch_cell_nums();
		__host__ void switch_phase();
	void init_with_file(const char* filename);
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
	float ca2p_avg_cache;
public:
	const CellIndex index;
	__device__ float& ageb(){ return ptr->ageb[index]; }
	__device__ float& agek(){ return ptr->agek[index]; }
	__device__ CellPos& pos(){ return ptr->pos[*ptr->current_phase][index]; }
	__device__ CellPos& next_pos(){ return ptr->pos[1 - *ptr->current_phase][index]; }
	__device__ CELL_STATE& state(){ return ptr->state[index]; }//__device__ CELL_STATE& state_ref(){ return ptr->state[index]; }
	__device__ CellIndex& fix_origin(){ return ptr->fix_origin[index]; }
	__device__ float& ca2p(){ return ptr->ca2p[*ptr->current_phase][index]; }
	__device__ float& next_ca2p(){ return ptr->ca2p[1-*ptr->current_phase][index]; }
	__device__ float ca2p_avg()const { return ca2p_avg_cache; }__device__ float& ca2p_avg_ref() { return ptr->ca2p_avg[index]; }
	__device__ float& ex_inert(){ return ptr->ex_inert[index]; }
	__device__ float& IP3(){ return ptr->IP3[*ptr->current_phase][index]; }
	__device__ float& next_IP3(){ return ptr->IP3[1-*ptr->current_phase][index]; }
	__device__ float& ex_fat(){ return ptr->ex_fat[index]; }
	__device__ float& in_fat(){ return ptr->in_fat[index]; }
	__device__ float& spr_nat_len(){ return ptr->spr_nat_len[index]; }
	__device__ int& rest_div_times(){ return ptr->rest_div_times[index]; }
	__device__ CellIndex& dermis_index(){ return ptr->dermis_index[index]; }
	__device__ CellDeviceWrapper(CellManager_Device* dev_ptr, CellIndex _index) :ptr(dev_ptr), index(_index){
		ca2p_avg_cache = dev_ptr->ca2p_avg[_index];
	}
	__device__ CellIndex& pair_index(){ return ptr->pair_index[index]; }
	__device__ CellConnectionData* connection_data_ptr(){ return &(ptr->connection_data[index]); }

	__device__ void remove(){
		ptr->remove_queue->push_back(index);
	}
};
float get_cell_zmax(CellManager*cm);
#endif /* CELLMANAGER_H_ */
