#include "define.h"
#include "CellManager.h"
#include "map.h"
#include "utils.h"
#define NSPHERE_COUNT (128)
#define MEMB_SPHERE_COUNT (32)
#define T_SPHERE_COUNT (768)

__constant__ char4 map_memb_sphere[MEMB_SPHERE_COUNT]; __constant__ int map_memb_count[1];
__constant__ char4 map_normal_sphere[NSPHERE_COUNT]; __constant__ int map_normal_count[1];
__constant__ char4 map_twice_sphere[T_SPHERE_COUNT]; __constant__ int map_twice_count[1];


void map_calc_init(){
	static const int _irx = (int)(R_max*NX / LX);
	static const int _iry = (int)(R_max*NY / LY);
	static const int _irz = (int)(R_max*NZ / LZ);
	char4 nm[1024];
	static const int m_irx = (int)(R_memb*NX / LX);
	static const int m_iry = (int)(R_memb*NY / LY);
	static const int m_irz = (int)(R_memb*NZ / LZ);

	static const int t_irx = (int)(2.0*R_max*NX / LX);
	static const int t_iry = (int)(2.0*R_max*NY / LY);
	static const int t_irz = (int)(2.0*R_max*NZ / LZ);

	int count = 0;
	for (char i = -_irx; i <= _irx; i++){
		float x = dx*i;
		for (char j = -_iry; j <= _iry; j++){
			float y = dy*j;
			for (char k = -_irz; k <= _irz; k++){
				float z = dz*k;
				if (x*x + y*y + z*z < R_max*R_max){
					nm[count].x = i;
					nm[count].y = j;
					nm[count].z = k;
					count++;
				}
			}
		}
	}
	cudaMemcpyToSymbol(map_normal_sphere, nm, sizeof(char4) * NSPHERE_COUNT,0,cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(map_normal_count, &count, sizeof(int), 0,cudaMemcpyHostToDevice);

	count = 0;
	for (char i = -m_irx; i <= m_irx; i++){
		float x = dx*i;
		for (char j = -m_iry; j <= m_iry; j++){
			float y = dy*j;
			for (char k = -m_irz; k <= m_irz; k++){
				float z = dz*k;
				if (x*x + y*y + z*z < R_memb*R_memb){
					nm[count].x = i;
					nm[count].y = j;
					nm[count].z = k;
					count++;
				}
			}
		}
	}

	cudaMemcpyToSymbol(map_memb_sphere, nm, sizeof(char4) * MEMB_SPHERE_COUNT, 0,cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(map_memb_count, &count, sizeof(int), 0,cudaMemcpyHostToDevice);

	count = 0;
	for (char i = -t_irx; i <= t_irx; i++){
		float x = dx*i;
		for (char j = -t_iry; j <= t_iry; j++){
			float y = dy*j;
			for (char k = -t_irz; k <= t_irz; k++){
				float z = dz*k;
				if (x*x + y*y + z*z < 4 * R_max*R_max){
					nm[count].x = i;
					nm[count].y = j;
					nm[count].z = k;
					count++;
				}
			}
		}
	}

	cudaMemcpyToSymbol(map_twice_sphere, nm, sizeof(char4) * T_SPHERE_COUNT, 0,cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(map_twice_count, &count, sizeof(int), 0,cudaMemcpyHostToDevice);

}

__global__ void setup_map_impl4(int ncell,int nmemb, CellPos* cs, unsigned int* cstate, int* cmap1_3d, float* cmap2_3d){
	const int i = blockIdx.x*blockDim.x + threadIdx.x;

	if (i < ncell){
		/*
		__shared__ char4 s_map_memb_sphere[MEMB_SPHERE_COUNT]; __shared__ int s_map_memb_count;
		__shared__ char4 s_map_normal_sphere[NSPHERE_COUNT]; __shared__ int s_map_normal_count;
		__shared__ char4 s_map_twice_sphere[T_SPHERE_COUNT]; __shared__ int s_map_twice_count;
		if (threadIdx.x == 0){
			s_map_memb_count = *map_memb_count;
			s_map_normal_count = *map_normal_count;
			s_map_twice_count = *map_twice_count;
			memcpy(s_map_memb_sphere, map_memb_sphere, sizeof(char4) * MEMB_SPHERE_COUNT);
			memcpy(s_map_normal_sphere, map_normal_sphere, sizeof(char4) * NSPHERE_COUNT);
			memcpy(s_map_twice_sphere, map_twice_sphere, sizeof(char4) * T_SPHERE_COUNT);
		}
		*/
		const unsigned int state = cstate[i];
		const CellPos cme = cs[i];
		
		const int4 lat = {
			(int)rintf(cme.x*inv_dx),
			(int)rintf(cme.y*inv_dy),
			(int)rintf(cme.z*inv_dz),
			0
		};
		__syncthreads();

		if (i<nmemb){
			for (int j = 0; j < map_memb_count[0]; j++){
				int ipx = (lat.x +map_memb_sphere[j].x + NX) % NX;
				int ipy = (lat.y + map_memb_sphere[j].y + NY) % NY;
				int ipz = lat.z + map_memb_sphere[j].z;

				cmap1_3d[midx<NX+1,NY+1,NZ+1>(ipx,ipy,ipz)] = ipz >= 0 && ipz < NZ?i:0x80808080;
			}
		}
		else{
			for (int j = 0; j < map_normal_count[0]; j++){
				int ipx = (lat.x + map_normal_sphere[j].x + NX) % NX;
				int ipy = (lat.y + map_normal_sphere[j].y + NY) % NY;
				int ipz = lat.z + map_normal_sphere[j].z;

				cmap1_3d[midx<NX + 1, NY + 1, NZ + 1>(ipx, ipy, ipz)] = ipz >= 0 && ipz < NZ ? i : 0x80808080;
			}
			const bool afm_cell = state == ALIVE || state == FIX || state == MUSUME;

			if (afm_cell){

				for (int j = 0; j < map_twice_count[0]; j++){
					int ipx = (lat.x + map_twice_sphere[j].x + NX) % NX;
					int ipy = (lat.y + map_twice_sphere[j].y + NY) % NY;
					int ipz = lat.z + map_twice_sphere[j].z;

					cmap2_3d[midx<NX + 1, NY + 1, NZ + 1>(ipx, ipy, ipz)] = ipz >= 0 && ipz < NZ ? 1.0f : 0.0f;
				}
			}
		}



	}
}
#define PAR_NUM 4

__global__ void setup_map_memb_impl(int nmemb, CellPos* cs, int* cmap1_3d, float* cmap2_3d){
	const int virt_id = threadIdx.x / PAR_NUM;
	const int par_offset = threadIdx.x%PAR_NUM;
	const int i = blockIdx.x*blockDim.x + virt_id;

	if (i < nmemb){
		/*
		__shared__ char4 s_map_memb_sphere[MEMB_SPHERE_COUNT]; __shared__ int s_map_memb_count;
		__shared__ char4 s_map_normal_sphere[NSPHERE_COUNT]; __shared__ int s_map_normal_count;
		__shared__ char4 s_map_twice_sphere[T_SPHERE_COUNT]; __shared__ int s_map_twice_count;
		if (threadIdx.x == 0){
		s_map_memb_count = *map_memb_count;
		s_map_normal_count = *map_normal_count;
		s_map_twice_count = *map_twice_count;
		memcpy(s_map_memb_sphere, map_memb_sphere, sizeof(char4) * MEMB_SPHERE_COUNT);
		memcpy(s_map_normal_sphere, map_normal_sphere, sizeof(char4) * NSPHERE_COUNT);
		memcpy(s_map_twice_sphere, map_twice_sphere, sizeof(char4) * T_SPHERE_COUNT);
		}
		*/
		//const unsigned int state = cstate[i];
		const CellPos cme = cs[i];

		const int3 lat = {
			(int)rintf(cme.x*inv_dx),
			(int)rintf(cme.y*inv_dy),
			(int)rintf(cme.z*inv_dz)
		};
		//__syncthreads();

		for (int j = par_offset; j < map_memb_count[0]; j += PAR_NUM){
				const int ipx = (lat.x + map_memb_sphere[j].x + NX) % NX;
				const int ipy = (lat.y + map_memb_sphere[j].y + NY) % NY;
				const int ipz = lat.z + map_memb_sphere[j].z;

				cmap1_3d[ipx*(NY + 1)*(NZ + 1) + ipy*(NZ + 1) + ipz] = ipz >= 0 && ipz < NZ ? i : 0x80808080;
		}




	}
}



__global__ void setup_map_non_memb_impl(int ncell, int offset, CellPos* cs, unsigned int* cstate, int* cmap1_3d, float* cmap2_3d){
	const int virt_id = threadIdx.x / PAR_NUM;
	const int par_offset = threadIdx.x%PAR_NUM;
	const int i = blockIdx.x*blockDim.x + virt_id + offset;

	if (i < ncell){
		const unsigned int state = cstate[i];
		const CellPos cme = cs[i];

		const int3 lat = {
			(int)rintf(cme.x*inv_dx),
			(int)rintf(cme.y*inv_dy),
			(int)rintf(cme.z*inv_dz)
		};
		//__syncthreads();

		
			for (int j = par_offset; j < map_normal_count[0]; j+=PAR_NUM){
				const int ipx = (lat.x + map_normal_sphere[j].x + NX) % NX;
				const int ipy = (lat.y + map_normal_sphere[j].y + NY) % NY;
				const int ipz = lat.z + map_normal_sphere[j].z;

				if (ipz >= 0 && ipz < NZ)cmap1_3d[midx<NX + 1, NY + 1, NZ + 1>(ipx, ipy, ipz)] = i;
			}
			const bool afm_cell = state == ALIVE || state == FIX || state == MUSUME;

			if (afm_cell){

				for (int j = par_offset; j < map_twice_count[0]; j+=PAR_NUM){
					const int ipx = (lat.x + map_twice_sphere[j].x + NX) % NX;
					const int ipy = (lat.y + map_twice_sphere[j].y + NY) % NY;
					const int ipz = lat.z + map_twice_sphere[j].z;

					if (ipz >= 0 && ipz < NZ)cmap2_3d[midx<NX + 1, NY + 1, NZ + 1>(ipx, ipy, ipz)] = 1.0f;
				}
			}
		



	}
}

void setup_map(CellManager* cm, int* cmap1_3d, float* cmap2_3d){
	cudaMemset(cmap1_3d, 0x80, sizeof(int)*(NX + 1)*(NY + 1)*(NZ + 1));
	cudaMemset(cmap2_3d, 0, sizeof(float)*(NX + 1)*(NY + 1)*(NZ + 1));

	setup_map_memb_impl << <cm->nmemb_host / 256 + 1, 256 * PAR_NUM >> >(cm->nmemb_host, cm->current_pos_host(), cmap1_3d, cmap2_3d);
	setup_map_non_memb_impl << <(cm->ncell_host - cm->nmemb_host) / 256 + 1, 256*PAR_NUM >> >(cm->ncell_host,cm->nmemb_host, cm->current_pos_host(),cm->state, cmap1_3d, cmap2_3d);
	cudaDeviceSynchronize();
}