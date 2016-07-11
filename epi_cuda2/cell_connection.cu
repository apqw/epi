/*
* cell_connection.cu
*
*  Created on: 2016/07/07
*      Author: yasu7890v
*/
#include "cell_connection.h"
#include "CellManager.h"
#include "utils.h"
#include <memory>
#include <cstdio>
#include <cassert>
#include <cfloat>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#define AREA_GRID_ORIGINAL 2.0f

#define AREA_GRID	(AREA_GRID_ORIGINAL + 1e-7f)
#define ANX ((int)(LX / AREA_GRID_ORIGINAL + 0.5f))
#define ANY ((int)(LY / AREA_GRID_ORIGINAL + 0.5f))
#define ANZ ((int)(LZ / AREA_GRID_ORIGINAL))
#define N3 64
#define N2 400

#define SEARCH_GRID_DIM 4

#define SEARCH_GRID_NUM (SEARCH_GRID_DIM*SEARCH_GRID_DIM*SEARCH_GRID_DIM)
#define SEARCH_PAR_NUM 8
#define SEARCH_THREAD_UPPER_LIMIT 64

#define GET_AREA_INFO_IDX(i) ((i)&0x007fffff)
#define GET_AREA_INFO_IDX_BY_PTR(ptr) GET_AREA_INFO_IDX(*reinterpret_cast<const unsigned int*>(ptr))
#define SET_AREA_INFO_MEMB_FLG(i) ((i)|1u<<23)
#define MERGE_AREA_INFO_STORE_COUNT(count,base) ((base)|(count<<24))
#define GET_AREA_INFO_STORE_COUNT(ptr) ((*reinterpret_cast<const unsigned int*>(ptr))>>24)
#define GET_AREA_INFO_MEMB_OR_NOT(i) ((i)&0x008000u!=0u)
#define GET_AREA_INFO_MEMB_OR_NOT_BY_PTR(ptr) GET_AREA_INFO_MEMB_OR_NOT(*reinterpret_cast<const unsigned int*>(ptr))


__global__ void init_grid(int ncell, CellPos* current_pos, CellManager_Device*const RESTRICT cmd, unsigned int*const RESTRICT aindx_3d, float4*const RESTRICT area_info){
	const unsigned int index = blockIdx.x*blockDim.x + threadIdx.x;

	if (index<ncell){
		//CellPos* current_pos=cmd->current_pos();

		const bool ismemb = cmd->state[index] == MEMB;
		__syncthreads();

		CellPos c = current_pos[index];
		//int aix, aiy, aiz;
		const int aix = ((int)rintf( (0.5f*LX - p_diff_x(0.5f*LX, c.x)) / AREA_GRID ))%ANX;
		const int aiy = ((int)rintf( (0.5f*LY - p_diff_y(0.5f*LY, c.y)) / AREA_GRID ))%ANY;
		const int aiz = ((int)rintf((min0(c.z)) / AREA_GRID))%ANZ;
		//printf("tesu: %d %d %d\n",aix,aiy,aiz);
		/*
		if ((aix >= ANX || aiy >= ANY || aiz >= ANZ || aix < 0 || aiy < 0 || aiz < 0)) {
			printf("err\n");
			printf("cx:%lf cy:%lf cz:%lf\n", c.x, c.y, c.z);
			printf("aix:%d aiy:%d aiz:%d\n", aix, aiy, aiz);
			assert(false);
		}
		*/

		//packing

		
		const int idx3d = midx<ANX, ANY, ANZ>(aix, aiy, aiz);
		const unsigned int listidx = atomicAdd(&aindx_3d[idx3d], 1);
		const int areaidx = idx3d*N3 + listidx;
		area_info[areaidx] = c;
		*(unsigned int*)(&area_info[areaidx].w) =
			ismemb ? SET_AREA_INFO_MEMB_FLG(index)
			: GET_AREA_INFO_IDX(index);

		//if (ismemb&&threadIdx.x == 25)printf("conn:%d\n", cmd->connection_data[index].connect_num);
		cmd->connection_data[index].connect_num = ismemb ? 8 : 0;

		if (listidx >= N3){
			printf("errroor:huge concentration num:%d aix:%d aiy:%d aiz:%d\n",listidx,aix,aiy,aiz);
		}

	}
}

__global__ void connect_proc(int nmemb, CellManager_Device*const RESTRICT cmd, const float4*const RESTRICT area_info){
	__shared__ int4 an;
	__shared__ CellPos my_pos;
	//__shared__ float my_radius;
	//radius always R_max
	__shared__ CellIndex my_index;
	if (threadIdx.x == 0){
		//printf("tesuya222\n");
		my_index = nmemb + blockIdx.x;
		my_pos = cmd->current_pos()[my_index];
		an.x = (((int)rintf(my_pos.x / AREA_GRID))+ANX)%ANX;
		an.y = (((int)rintf(my_pos.y / AREA_GRID))+ANY)%ANY;
		an.z = ((int)rintf(min0(my_pos.z) / AREA_GRID))%ANZ;
	}

	const int virtual_id = threadIdx.x / SEARCH_PAR_NUM;
	const int offset = threadIdx.x%SEARCH_PAR_NUM;
	int cz = virtual_id;
	//MEMO:���̂�SEARCH_PAR_NUM�Ŋ����ăn�}���Ă�
	int cy = ((int)(cz / SEARCH_GRID_DIM)); cz = cz%SEARCH_GRID_DIM;
	int cx = ((int)(cy / SEARCH_GRID_DIM)); cy = cy%SEARCH_GRID_DIM;
	cx = cx%SEARCH_GRID_DIM;
	__syncthreads();
	if (virtual_id <= SEARCH_GRID_NUM){
		const int aix = (an.x - (int)(SEARCH_GRID_DIM / 2) + cx + ANX) % ANX;
		const int aiy = (an.y - (int)(SEARCH_GRID_DIM / 2) + cy + ANY) % ANY;
		const int aiz = (an.z - (int)(SEARCH_GRID_DIM / 2) + cz + ANZ) % ANZ;
		const int idx3d = midx<ANX, ANY, ANZ>(aix, aiy, aiz)*N3;
		const int store_count = GET_AREA_INFO_STORE_COUNT(&area_info[idx3d + 0].w); //*(unsigned int*)&area_info[idx3d + 0].w>>24;
		
		for (int i = offset; i<store_count; i += SEARCH_PAR_NUM){
			const float4 cpos = area_info[idx3d + i];
			const CellIndex cidx = GET_AREA_INFO_IDX_BY_PTR(&cpos.w);
			const bool op_memb = GET_AREA_INFO_MEMB_OR_NOT_BY_PTR(&cpos.w);
			if (my_index>cidx){
				const float distSq = p_cell_dist_sq(my_pos, cpos);
				const float rad_sum = R_max + (op_memb ? R_memb : R_max);
				if (distSq<LJ_THRESH*LJ_THRESH*rad_sum*rad_sum){
					//if(cidx<1000)printf("cidx:%d\n", cidx);
					cmd->connection_data[my_index].add_atomic(cidx);
					cmd->connection_data[cidx].add_atomic(my_index);
				}
			}
		}
	}

}

__global__ void complete_area_info(const unsigned int*const RESTRICT aindx_3d, float4*const RESTRICT area_info){
	const int idx = midx<ANX, ANY, ANZ>(blockIdx.x, blockIdx.y, threadIdx.x);
	const unsigned int w = *reinterpret_cast<unsigned int*>(&area_info[idx*N3 + 0].w);
	const unsigned int num = aindx_3d[idx];
	*reinterpret_cast<unsigned int*>(&area_info[idx*N3 + 0].w)
		= MERGE_AREA_INFO_STORE_COUNT(num, 0x00ffffff&w);
}

__global__ void set_dermis(int ncell, int offset, CellManager_Device* cmd){
	const CellIndex index = blockIdx.x*blockDim.x + threadIdx.x + offset;
	if (index<ncell){
		const CELL_STATE my_state = cmd->state[index];
		if (my_state == FIX || my_state == MUSUME){
			const CellPos cme = cmd->current_pos()[index];
			const unsigned int cnum = cmd->connection_data[index].connect_num;
			int dermis_index = -1;
			float distSq = FLT_MAX;
			for (int i = 0; i<cnum; i++){
				const int opidx = cmd->connection_data[index].connect_index[i];
				const CellPos opcs = cmd->current_pos()[opidx];
				const CELL_STATE op_state = cmd->state[opidx];
				if (op_state == MEMB){
					const float tmp_dist_sq = p_cell_dist_sq(cme, opcs);
					if (tmp_dist_sq<distSq){
						dermis_index = opidx;
						distSq = tmp_dist_sq;
					}
				}
			}
			cmd->dermis_index[index] = dermis_index;
		}
	}
}


void dbg_conn_info(unsigned int* aindx,int count){
	static std::vector<unsigned int> tmp(ANX*ANY*ANZ);
	cudaMemcpy(&tmp[0], aindx, sizeof(unsigned int)*ANX*ANY*ANZ, cudaMemcpyDeviceToHost);
	std::ofstream ofs(("dbg_ci" + int_to_string(count)).c_str());
	for (int i = 0; i < ANX; i++){
		for (int j = 0; j < ANY; j++){
			for (int k = 0; k < ANZ; k++){
				ofs << i << "," << j << "," << k << " " << tmp[i*ANY*ANZ + j*ANZ + k] << std::endl;
			}
		}
	}

}
void connect_cell(CellManager* cm){
	static device_alloc_ctor<unsigned int> aindx(ANX*ANY*ANZ);
	static device_alloc_ctor<float4> area_info(ANX*ANY*ANZ*N3);
	static int count = 0;
	aindx.set_zero();
	area_info.set_zero();
	init_grid << <256, 256 >> >(cm->ncell_host, cm->current_pos_host(), cm->dev_ptr, aindx.ptr, area_info.ptr);
	cudaDeviceSynchronize();
	complete_area_info << <dim3(ANX, ANY), ANZ >> >(aindx.ptr, area_info.ptr);
	cudaDeviceSynchronize();
	connect_proc << <cm->ncell_host - cm->nmemb_host, SEARCH_THREAD_UPPER_LIMIT*SEARCH_PAR_NUM >> >(cm->nmemb_host, cm->dev_ptr, area_info.ptr);
	cudaDeviceSynchronize();
	set_dermis << <(cm->ncell_host - cm->nmemb_host - cm->nder_host) / 256 + 1, 256 >> >(cm->ncell_host, cm->nder_host + cm->nmemb_host, cm->dev_ptr);
	cudaDeviceSynchronize();
	cm->no_need_reconnect();
}
