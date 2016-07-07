#pragma once

#include <cuda_runtime.h>
#include "device_launch_parameters.h"
#include "helper_cuda.h"
using cell_pos_set = float4;
#define cdefd static const float
#define cdefi static const int
#define cdefui static const unsigned int
#define cdefs static const size_t

#define defd static double
#define rdefd static double&
#define defi static int
#define rdefi static int&
#define defui static unsigned int
#define rdefui static unsigned int&
#define defs static size_t
#define rdefs static size_t&

#define cmax(a,b) ((a)>(b)?(a):(b))
#ifdef _WIN32
#define RESTRICT __restrict
#else
#define RESTRICT __restrict__
#endif
#define MAX_CONNECT_CELL_NUM (200u)
#define NMX (300u)
#define NMY (150u)
#define MEMB_NUM (NMX*NMY)
#define MAX_NON_MEMB_NUM (20000u)
#define MAX_CELL_NUM (65536u)

#define LX (100.0f)
#define LY (50.0f)
#define LZ (100.0f)

#define COMPRESS_FACTOR (6u)

#define R_max (1.4f)
#define R_der (1.4f)
#define R_memb (1.0f)

#define NX (200)
#define NY (100)
#define NZ (200)

#define dx (LX/NX)
#define dy (LY/NY)
#define dz (LZ/NZ)

#define inv_dx (NX/LX)
#define inv_dy (NY/LY)
#define inv_dz (NZ/LZ)

enum CELL_STATE :unsigned int {
	ALIVE = 0,
	DEAD = 1,
	DISA = 2,
	UNUSED = 3, //
	FIX = 4,
	BLANK = 5,
	DER = 6,
	MUSUME = 7,
	AIR = 8,
	MEMB = 9
};

//warn:huge
struct connected_index_set{
	int connected_num=0;
	int index[MAX_CONNECT_CELL_NUM];
	float gj[MAX_CONNECT_CELL_NUM];
	//cell_pos_set cp[MAX_CONNECT_CELL_NUM];
};

//2 means new/old (this switches into old/new)
struct DeviceData{
	CELL_STATE* c_state_d;
	cell_pos_set* c_pos_d[2]; //w:=state(duplicated)
	int* c_fix_origin_d;
	connected_index_set* c_connected_index_d;
	int* c_pair_index_d;
	float* c_ca2p_d[2];
	float* c_ca2p_avg_d;
	float* c_ex_inert_d;
	float* c_IP3_d;
	float* c_agek_d;
	float* c_ageb_d;
	float* c_ex_fat_d;
	float* c_in_fat_d;
	float* c_spr_nat_len_d;
	//float radius[];
	int* c_rest_div_times_d;
	int* c_dermis_index_d;
	float* zzmax;
	int ncell = 0;
	int nder = 0;
	int nmemb = 0;
	int nother=0;
	int current=0;
	int* block_max_store;
	int* sw;

	DeviceData(){
		
		checkCudaErrors(cudaMalloc((void**)&c_state_d, sizeof(CELL_STATE)*MAX_CELL_NUM));
		checkCudaErrors(cudaMalloc((void**)&c_pos_d[0], sizeof(cell_pos_set)*MAX_CELL_NUM));
		checkCudaErrors(cudaMalloc((void**)&c_pos_d[1], sizeof(cell_pos_set)*MAX_CELL_NUM));
		checkCudaErrors(cudaMalloc((void**)&c_fix_origin_d, sizeof(int)*MAX_CELL_NUM));
		checkCudaErrors(cudaMalloc((void**)&c_connected_index_d, sizeof(connected_index_set)*MAX_CELL_NUM));
		checkCudaErrors(cudaMalloc((void**)&c_pair_index_d, sizeof(int)*MAX_CELL_NUM));
		checkCudaErrors(cudaMalloc((void**)&c_ca2p_d[0], sizeof(float)*MAX_CELL_NUM));
		checkCudaErrors(cudaMalloc((void**)&c_ca2p_d[1], sizeof(float)*MAX_CELL_NUM));
		checkCudaErrors(cudaMalloc((void**)&c_ca2p_avg_d, sizeof(float)*MAX_CELL_NUM));
		checkCudaErrors(cudaMalloc((void**)&c_ex_inert_d, sizeof(float)*MAX_CELL_NUM));
		checkCudaErrors(cudaMalloc((void**)&c_IP3_d, sizeof(float)*MAX_CELL_NUM));
		checkCudaErrors(cudaMalloc((void**)&c_agek_d, sizeof(float)*MAX_CELL_NUM));
		checkCudaErrors(cudaMalloc((void**)&c_ageb_d, sizeof(float)*MAX_CELL_NUM));
		checkCudaErrors(cudaMalloc((void**)&c_ex_fat_d, sizeof(float)*MAX_CELL_NUM));
		checkCudaErrors(cudaMalloc((void**)&c_in_fat_d, sizeof(float)*MAX_CELL_NUM));
		checkCudaErrors(cudaMalloc((void**)&c_spr_nat_len_d, sizeof(float)*MAX_CELL_NUM));
		checkCudaErrors(cudaMalloc((void**)&c_rest_div_times_d, sizeof(int)*MAX_CELL_NUM));
		checkCudaErrors(cudaMalloc((void**)&c_dermis_index_d, sizeof(int)*MAX_CELL_NUM));
		checkCudaErrors(cudaMalloc((void**)&zzmax, sizeof(float)));
		checkCudaErrors(cudaMalloc((void**)&block_max_store, sizeof(int)*512));
		checkCudaErrors(cudaMalloc((void**)&sw, sizeof(int)));//plz set!!
	}
	~DeviceData(){
		cudaFree(c_state_d);
		cudaFree(c_pos_d[0]);
		cudaFree(c_pos_d[1]);
		cudaFree(c_fix_origin_d);
		cudaFree(c_connected_index_d);
		cudaFree(c_pair_index_d);
		cudaFree(c_ca2p_d[0]);
		cudaFree(c_ca2p_d[1]);
		cudaFree(c_ca2p_avg_d);
		cudaFree(c_ex_inert_d);
		cudaFree(c_IP3_d);
		cudaFree(c_agek_d);
		cudaFree(c_ageb_d);
		cudaFree(c_ex_fat_d);
		cudaFree(c_in_fat_d);
		cudaFree(c_spr_nat_len_d);
		cudaFree(c_rest_div_times_d);
		cudaFree(c_dermis_index_d);
		cudaFree(zzmax);
		cudaFree(sw);
		cudaFree(block_max_store);
	}
};

#define SYSTEM 0
#define WHOLE 0
#define BASAL 1

namespace cont {


	cdefs STATE_NUM = 10;
	__constant__ cdefd DT_Cell = 0.01;
	__constant__ cdefd DT_Ca = 0.01;//0.02;
	//cdefi MALIG_NUM = 0;

	/*
	LX,LY,LZ:ŒvŽZ—Ìˆæ‚ÌƒTƒCƒY
	*/
	//__constant__ cdefd LX = 100.0;
	//__constant__ cdefd LY = 50.0;
	//__constant__ cdefd LZ = 100.0;

	/*
	NX,NY,NZ:ŒvŽZ—Ìˆæ‚Ì•ªŠ„”
	*/

	/*
	dx,dy,dz:ŒvŽZ—Ìˆæ‚Ì•ªŠ„•(L**,N**‚©‚çŒvŽZ‚³‚ê‚é)
	*/


	
	//cdefd COMPRESS_FACTOR = 6;//ok
	cdefui MEMB_ADHE_RANGE = (unsigned int)cmax((int)((R_max / R_der)*COMPRESS_FACTOR / 2) - 2, 0);//test

	__constant__ cdefd THRESH_SP = 3.0; //ok 3.0->1.0

	__constant__ cdefd LJ_THRESH = 1.2;//ok
	cdefui MEMB_ADJ_CONN_NUM = 8;



	__constant__ cdefd THRESH_DEAD = 22.0;//ok

	__constant__ cdefui NUM_ITR = 4 * ((int)1e6); //twice



	__constant__ cdefd Ca_avg_time = 10.0;


	__constant__ cdefd ca2p_init = 0.122;

	__constant__ cdefd ex_inert_init = 0.97;

	__constant__ cdefd IP3_init = 0;

	__constant__ cdefd gj_init = 0.99;

	__constant__ cdefd ATP_init = 0;

	__constant__ cdefd ext_stim_init = 0;



	cdefd FAC_MAP = 2.0;//ok

	cdefi SW_THRESH = 20;//ok

	cdefd T_TURNOVER = 6000.0;//
	cdefi NUM_SC_INIT = 1;//ok ha/???->19->1
	



}

enum DIRECTION :unsigned int{
	U = 0, L = 1, B = 2, R = 3
};//why not B->D?

#define STOCHASTIC (1)
#define M_PI_F (3.141592654f)
#define MALIGNANT (0)
#define stoch_corr_coef (1.0f)