#pragma once
#undef NDEBUG
#include "global_cuda.h"
#include <device_functions.h>
#include "cell_connection_cuda.h"
#include <stdio.h>
#include "utils_cuda.cuh"
#include <cassert>
//#include <thrust/sort.h>
#define AREA_GRID_ORIGINAL 2.0//ok

/** グリッドサイズ (?)*/
#define AREA_GRID (AREA_GRID_ORIGINAL + 1e-7)//ok



/** X方向のグリッド数*/
__constant__  static const int ANX = (LX / AREA_GRID_ORIGINAL + 0.5);//ok
/** Y方向のグリッド数*/
__constant__ static const  int ANY = (LY / AREA_GRID_ORIGINAL + 0.5);//ok
/** Z方向のグリッド数*/
__constant__ static const  int ANZ = (LZ / AREA_GRID_ORIGINAL);//ok
/** グリッド1つ当たりの細胞格納数上限 */
__constant__ static const  int N3 = 64; //max grid cell num //ok
/** 1細胞の接続最大数 */
__constant__ static const  int N2 = 400; //max conn num //ok
#define SEARCH_GRID_DIM (5)
#define EXTRACT_STATE_CS(cs) (*(unsigned int*)(&(cs).w))

__global__ void init_grid(int ncell, unsigned int* cstate, cell_pos_set* cs, connected_index_set*cis, unsigned int* aindx_3d, float4* area_info){
	unsigned int index = blockIdx.x*blockDim.x + threadIdx.x;
	//while (true){ index++; }
	//exit(1);
	if (index < ncell){
		cell_pos_set c = cs[index];
		bool ismemb = cstate[index] == MEMB;
		using namespace cont;
		int aix = (int)((0.5f*LX - p_diff_x(0.5f*LX, c.x)) / (float)AREA_GRID);
		int aiy = (int)((0.5f*LY - p_diff_y(0.5f*LY, c.y)) / (float)AREA_GRID);
		int aiz = (int)((min0(c.z)) / (float)AREA_GRID);


		if ((aix >= ANX || aiy >= ANY || aiz >= ANZ || aix < 0 || aiy < 0 || aiz < 0)) {
			printf("err\n");
			printf("cx:%lf cy:%lf cz:%lf\n", c.x, c.y, c.z);
			printf("aix:%d aiy:%d aiz:%d\n", aix, aiy, aiz);
			assert(false);
		}

		/*
		area_info:
		x=x pos
		y=y pos
		z=z pos
		w=
		head 8bit =connection num(unsigned char)
		next 1bit = is memb
		rest 23bit = index(as unsigned32bit int)



		*/

		int idx3d = aix*ANY*ANZ + aiy*ANZ + aiz;
		unsigned int listidx = atomicAdd(&aindx_3d[idx3d], 1);
		int areaidx = idx3d*N3 + listidx;
		area_info[areaidx] = c;
		unsigned int wvalue = ismemb ? index | 1u << 23 : index&(0x007fffff);

		*(unsigned int*)&area_info[areaidx].w = wvalue;
		cis[index].connected_num = ismemb ? 8 : 0;
		/*
		このassertは消してよい
		*/
		if (listidx >= N3){
			printf("error:aaaaaa\n");
		}
		//assert(aindx[aix][aiy][aiz] < (cint)N3);
		//c->connected_cell.force_set_count(c->state == MEMB ? MEMB_ADJ_CONN_NUM : 0);
	}
}
#define ADD_CELL_CONNECT(cis,me,add) cis[me].index[atomicAdd(&cis[me].connected_num,1)]=add
__global__ void connect_proc(int n_cell, int n_memb, cell_pos_set* cs, connected_index_set* cis, float4* area_info){

	using namespace cont;
	__shared__ int4 an;
	__shared__ cell_pos_set my_pos;
	__shared__ float my_radius;
	//__shared__ int conn_num[128];
	//__shared__ int conn_count;
	int my_index = n_memb + blockIdx.x;

	if (threadIdx.x == 0){
		/*
		for (int jj = 0; jj < 128; jj++){
			conn_num[jj] = -1;
		}
		*/
		//conn_count = 0;
		my_pos = cs[my_index];
		my_radius = get_radius(EXTRACT_STATE_CS(my_pos));
		an.x = my_pos.x / (float)AREA_GRID;
		an.y = my_pos.y / (float)AREA_GRID;
		an.z = my_pos.z / (float)AREA_GRID;
		//if(my_index>=n_cell-1)
	}
#define DIV_NUM 8
	int virt_id = threadIdx.x / DIV_NUM;
	int divflg = threadIdx.x % DIV_NUM;
	int cz = virt_id;
	int cy = ((int)(cz / (int)SEARCH_GRID_DIM)); cz = cz%SEARCH_GRID_DIM;
	int cx = ((int)(cy / (int)SEARCH_GRID_DIM)) % SEARCH_GRID_DIM; cy = cy%SEARCH_GRID_DIM;
	__syncthreads();//do work before sync as much as possible

	if (virt_id <= SEARCH_GRID_DIM*SEARCH_GRID_DIM*SEARCH_GRID_DIM){

		int aix = an.x - (SEARCH_GRID_DIM / 2) + cx; aix = (aix + ANX) % ANX;
		int aiy = an.y - (SEARCH_GRID_DIM / 2) + cy; aiy = (aiy + ANY) % ANY;
		int aiz = an.z - (SEARCH_GRID_DIM / 2) + cz; aiz = (aiz + ANZ) % ANZ;
		//printf("aix:%d aiy:%d aiz:%d\n", aix, aiy, aiz);
		int idx3d = (aix*ANY*ANZ + aiy*ANZ + aiz)*N3;
		int m = *(unsigned int*)&area_info[idx3d + 0].w >> 24;

		for (int i = divflg; i < m; i += DIV_NUM){
			float4 cpos = area_info[idx3d + i];
			//__syncthreads();
			int cidx = (*(unsigned int*)&cpos.w)&(0x007fffff);
			//if (threadIdx.x == 25)printf("idx:%d\n",cidx);
			bool op_memb = (*(unsigned int*)&cpos.w)&(0x00800000) != 0;
			if (my_index > cidx){
				float diffx = p_diff_x(my_pos.x, cpos.x);
				float diffy = p_diff_y(my_pos.y, cpos.y);
				float diffz = my_pos.z - cpos.z;
				float rad_sum = my_radius + (op_memb ? R_memb : R_max);
				if (diffx*diffx + diffy*diffy + diffz*diffz < LJ_THRESH*LJ_THRESH*rad_sum*rad_sum){

					ADD_CELL_CONNECT(cis, my_index, cidx);
					ADD_CELL_CONNECT(cis, cidx, my_index);
				}
			}
		}

	}
}

__global__ void connect_proc2(int n_cell, int n_memb, cell_pos_set* cs, unsigned int* cstate, connected_index_set* cis, float4* area_info){

	using namespace cont;
	int index = blockIdx.x*blockDim.x + threadIdx.x;
	__shared__ int conn_num[256];
	if (index < n_cell){
		int mylist[128];
		int mylist_count = 0;
		cell_pos_set my_pos = cs[index];
		bool me_memb = cstate[index] == MEMB;
		__syncthreads();
		float my_radius = EXTRACT_STATE_CS(my_pos);
		int4 an = { (int)(my_pos.x / (float)AREA_GRID),
			(int)(my_pos.y / (float)AREA_GRID),
			(int)(my_pos.z / (float)AREA_GRID), 0 };

		for (int i = an.x - 2; i <= an.x + 2; i++){
			int aix = (i + ANX) % ANX;
			for (int j = an.y - 2; j <= an.y + 2; j++){
				int aiy = (j + ANY) % ANY;
				for (int k = an.z - 2; k <= an.z + 2; k++){
					int aiz = (k + ANZ) % ANZ;
					int idx3d = (aix*ANY*ANZ + aiy*ANZ + aiz)*N3;
					int m = *(unsigned int*)&area_info[idx3d + 0].w >> 24;
					for (int l = 0; l < m; l++){
						float4 cpos = area_info[idx3d + l];
						int cidx = (*(unsigned int*)&cpos.w)&(0x007fffff);
						bool op_memb = (*(unsigned int*)&cpos.w)&(0x00800000) != 0;
						if (index != cidx && !(op_memb&&me_memb)){
							float diffx = p_diff_x(my_pos.x, cpos.x);
							float diffy = p_diff_y(my_pos.y, cpos.y);
							float diffz = my_pos.z - cpos.z;
							float rad_sum = my_radius + (op_memb ? R_memb : R_max);
							if (diffx*diffx + diffy*diffy + diffz*diffz < LJ_THRESH*LJ_THRESH*rad_sum*rad_sum){
								//ADD_CELL_CONNECT(cis, index, cidx);
								//mylist[mylist_count++] = cidx;
								//ADD_CELL_CONNECT(cis, cidx, index);
								mylist[mylist_count] = cidx;
								mylist_count++;
							}
						}

					}
				}
			}
		}

		cis[index].connected_num = mylist_count;
		/*
		for (int kk = 0; kk < mylist_count; kk++){
		cis[index].index[kk] = mylist[kk];
		}
		*/
	}
}

__global__ void complete_area_info(unsigned int* aindx_3d, float4* area_info){
	int idx = blockIdx.x*ANY*ANZ + blockIdx.y*ANZ + threadIdx.x;
	unsigned int w = *(unsigned int*)&area_info[idx*N3 + 0].w;
	unsigned int num = aindx_3d[idx];

	*(unsigned int*)&area_info[idx*N3 + 0].w = ((unsigned int)(0x00ffffff) & w) | (num << 24);
}

__global__ void set_dermis(int offset, int ncell, cell_pos_set* pos, connected_index_set* cis, int* dermis_arr){
	int index = blockIdx.x*blockDim.x + threadIdx.x + offset;
	if (index < ncell){
		cell_pos_set cme = pos[index];
		int cnum = cis[index].connected_num;
		int dermis_index = -1;
		float distSq = FLT_MAX;
		float tmpdistSq = 0;
		for (int i = 0; i < cnum; i++){
			int opidx = cis[index].index[i];
			cell_pos_set opcs = pos[cis[index].index[i]];
			if (*(unsigned int*)&opcs.w == MEMB){
				tmpdistSq = p_cell_dist_sq_d(cme, opcs);
				if (tmpdistSq < distSq){
					dermis_index = opidx;
				}
			}
		}
		dermis_arr[index] = dermis_index;
	}
}

__global__ void conn_sort(int ncell, int nmemb, unsigned int* cstate, connected_index_set* cis){
	int index = blockIdx.x*blockDim.x + threadIdx.x;
	if (index < ncell){
		int offset = cstate[index] == MEMB ? 8 : 0;
		int bucket[10][128] = {}; int bucket_count[10] = { 0 };
		for (int i = offset; i < cis[index].connected_num; i++) {
			int idx = cis[index].index[i];
			unsigned int state = cstate[idx];
			//bucket[state][bucket_count[state]++] = idx;
		}
		/*
		int count = 0;
		for (int i = 0; i < 10; i++){
			for (int j = 0; j < bucket_count[i]; j++){
				cis[index].index[offset + (count++)] = bucket[i][j];
			}
		}
		*/
	}
}

void connect_cell(DeviceData* d){

	struct uint_alloc{
		unsigned int* ptr;
		uint_alloc(size_t size){
			checkCudaErrors(cudaMalloc((void**)&ptr, size));
		}
		~uint_alloc(){
			cudaFree(ptr);
		}
	};

	struct f4_alloc{
		float4* ptr;
		f4_alloc(size_t size){
			checkCudaErrors(cudaMalloc((void**)&ptr, size));
		}
		~f4_alloc(){
			cudaFree(ptr);
		}
	};
	static uint_alloc aindx(sizeof(unsigned int)*ANX*ANY*ANZ);
	//static uint_alloc area(sizeof(unsigned int)*ANX*ANY*ANZ*N3);
	static f4_alloc area_info(sizeof(float4)*ANX*ANY*ANZ*N3);
	cudaMemset(aindx.ptr, 0, sizeof(unsigned int)*ANX*ANY*ANZ);
	cudaMemset(area_info.ptr, 0, sizeof(float4)*ANX*ANY*ANZ*N3);
	init_grid << <256, 256 >> >(d->ncell, (unsigned int*)d->c_state_d, d->c_pos_d[d->current], d->c_connected_index_d, aindx.ptr, area_info.ptr);
	cudaDeviceSynchronize();
	dim3 grd(ANX, ANY);

	complete_area_info << <grd, ANZ >> >(aindx.ptr, area_info.ptr);
	cudaDeviceSynchronize();
	connect_proc << <d->ncell - d->nmemb, 128*8 >> >(d->ncell, d->nmemb, d->c_pos_d[d->current], d->c_connected_index_d, area_info.ptr);
	cudaDeviceSynchronize();
	//conn_sort << < d->ncell / 256 + 1, 256 >>> (d->ncell, d->nmemb, (unsigned int*)d->c_state_d, d->c_connected_index_d);
	//cudaDeviceSynchronize();
	set_dermis << <(d->ncell - d->nmemb - d->nder) / 256 + 1, 256 >> >(d->nmemb + d->nder, d->ncell, d->c_pos_d[d->current], d->c_connected_index_d, d->c_dermis_index_d);
	cudaDeviceSynchronize();
}