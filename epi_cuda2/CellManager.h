/*
 * CellManager.h
 *
 *  Created on: 2016/07/07
 *      Author: yasu7890v
 */

#ifndef CELLMANAGER_H_
#define CELLMANAGER_H_
#include "define.h"

struct CellConnectionData{
	unsigned int connect_num;
	CellIndex connect_index[MAX_CONNECT_CELL_NUM];
	float gj[MAX_CONNECT_CELL_NUM];
	__device__ void add_atomic(CellIndex idx);
};

struct CellDataStruct{

	int* __block_max_store;
	int* current_phase;
	//__global__ void set_connection_data_ref

	CELL_STATE* state;
	CellPos* pos[2];
	CellIndex* fix_origin;
	CellConnectionData* connection_data;
	CellIndex* pair_index;
	FloatArr* ca2p[2];
	FloatArr ca2p_avg;
	FloatArr ex_inert;
	FloatArr IP3;
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
};

class CellManager_Device:public CellDataStruct{
public:
	__device__ CellPos* current_pos(){
		return pos[*(current_phase)];
	}
};

class CellManager:public CellDataStruct{

	void alloc();
	void dealloc();
	void memb_init(int,CellConnectionData[]);
public:
	int current_phase_host;
	int ncell_host;
	int nder_host;
	int nmemb_host;
	__host__ void set_cell_nums(int _ncell,int _nmemb,int _nder);
		//__host__ void fetch_cell_nums();
		__host__ void switch_phase();
	void init_with_file(const char* filename);
	CellPos* current_pos_host();
	CellManager_Device* dev_ptr;
	CellManager():current_phase_host(0),ncell_host(0),nder_host(0),nmemb_host(0){
		alloc();
	}
	~CellManager(){
		dealloc();
	}
	static const int DEDUCT_ZMAX_DIV_NUM=8;
};

#endif /* CELLMANAGER_H_ */
