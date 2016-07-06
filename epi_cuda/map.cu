#include "global_cuda.h"
#include "utils_cuda.cuh"
#include "map.h"


#define NSPHERE_COUNT (128)
#define MEMB_SPHERE_COUNT (32)
#define T_SPHERE_COUNT (768)
__constant__ static const int irx = (int)(R_max*NX / LX);
__constant__ static const int iry = (int)(R_max*NY / LY);
__constant__ static const int irz = (int)(R_max*NZ / LZ);
__constant__ char4 map_memb_sphere[MEMB_SPHERE_COUNT]; __constant__ int map_memb_count;
__constant__ char4 map_normal_sphere[NSPHERE_COUNT]; __constant__ int map_normal_count;
__constant__ char4 map_twice_sphere[T_SPHERE_COUNT]; __constant__ int map_twice_count;

//need init
__global__ void setup_map_impl(int ncell, cell_pos_set* cs, unsigned int* cstate, int* cmap1_3d, float* cmap2_3d){
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i < ncell){
		unsigned int state = cstate[i];
		cell_pos_set cme = cs[i];
		const bool afm_cell = state == ALIVE || state == FIX || state == MUSUME;
		int4 lat = {
			((int)(cme.x*inv_dx + NX)) % NX,
			((int)(cme.y*inv_dy + NY)) % NY,
			((int)(cme.z*inv_dz + NZ)) % NZ,
			0
		};


		if (!afm_cell){
			const int z_b = lat.z - irz;
			const int z_u = lat.z + irz;
			const int zmin = z_b >= 1 ? z_b : 1;
			const int zmax = z_u < (int)NZ ? z_u : (int)NZ - 1;
			float rad = get_radius(state);
			for (int k = lat.x - irx; k <= lat.x + irx; k++){
				float mx = k*dx;
				int ipx = (k+NX)%NX;
				float diffx = p_diff_x(mx, cme.x);
				float diffxSq = diffx*diffx;

				for (int l = lat.y - iry; l <= lat.y + iry; l++){
					float my = l*dy;
					int ipy = (l + NY) % NY;
					float diffy = p_diff_y(my, cme.y);
					float diffySq = diffy*diffy;
					for (int ipz = zmin; ipz <= zmax; ipz++){
						float mz = ipz*dz;
						float diffz = mz - cme.z;
						if (diffz*diffz + diffySq + diffxSq < rad*rad){
							cmap1_3d[ipx*(NY + 1)*(NZ + 1) + ipy*(NZ + 1) + ipz] = i;
						}
					}
				}
			}
		}
		
		else{
			int _irz = 2 * irz;
			int _iry = 2 * iry;
			int _irx = 2 * irx;
			const int z_b = lat.z - _irz;
			const int z_u = lat.z + _irz;
			const int zmin = z_b >= 1 ? z_b : 1;
			const int zmax = z_u < (int)NZ ? z_u : (int)NZ - 1;
			float rad = get_radius(state);
			for (int k = lat.x - _irx; k <= lat.x + _irx; k++){
				float mx = k*dx;
				int ipx = (k + NX) % NX;
				float diffx = p_diff_x(mx, cme.x);
				float diffxSq = diffx*diffx;

				for (int l = lat.y - _iry; l <= lat.y + _iry; l++){
					float my = l*dy;
					int ipy = (l + NY) % NY;
					float diffy = p_diff_y(my, cme.y);
					float diffySq = diffy*diffy;
					for (int ipz = zmin; ipz <= zmax; ipz++){
						float mz = ipz*dz;
						float diffz = mz - cme.z;
						float distSq = diffz*diffz + diffySq + diffxSq;
						if (distSq < 4*rad*rad){
							cmap2_3d[ipx*(NY + 1)*(NZ + 1) + ipy*(NZ + 1) + ipz] = 1.0f;
							if (distSq < rad*rad){
								cmap1_3d[ipx*(NY + 1)*(NZ + 1) + ipy*(NZ + 1) + ipz] = i;
							}
						}
					}
				}
			}
		}
		
	}
}

__global__ void setup_map_impl2(int ncell, cell_pos_set* cs, unsigned int* cstate, int* cmap1_3d, float* cmap2_3d){
	int i = blockIdx.x;
	__shared__ unsigned int state;
	__shared__ cell_pos_set cme;
	__shared__ bool afm_cell;
	__shared__ int4 lat;
	__shared__ int xmin;
	__shared__ int xmax;
	__shared__ int ymin;
	__shared__ int ymax;
	__shared__ int zmin;
	__shared__ int zmax;
	if (threadIdx.x == 0){
		state = cstate[i];
		cme = cs[i];
		afm_cell = state == ALIVE || state == FIX || state == MUSUME;
		lat.x = ((int)(cme.x*inv_dx + NX)) % NX;
		lat.y = ((int)(cme.y*inv_dy + NY)) % NY;
		lat.z = ((int)(cme.z*inv_dz + NZ)) % NZ;

		int z_b = lat.z - irz;
		int z_u = lat.z + irz;
		zmin = z_b >= 1 ? z_b : 1;
		zmax = z_u < (int)NZ ? z_u : (int)NZ - 1;
		xmin = lat.x - irx; xmax = lat.x + irx;
		ymin = lat.y - iry; ymax = lat.y + iry;

	}
	__syncthreads();
	int cz = threadIdx.x;
	int cy = ((int)(cz / (int)5)); cz = cz%(zmax-zmin);
	int cx = ((int)(cy / (int)5)) % 5; cy = cy%5;
	

	if (threadIdx.x <= 25 * (zmax - zmin)){
		int ipx = (xmin + cx+NX)%NX;
		int ipy = (ymin + cy+NY)%NY;
		int ipz = zmin + cz;
		float mx = ipx*dx;
		float my = ipy*dy;
		float mz = ipz*dz;
		float diffx = p_diff_x(mx, cme.x);
		float diffy = p_diff_y(my, cme.y);
		float diffz = mz - cme.z;

		if (diffx*diffx + diffy*diffy + diffz*diffz < 1.4*1.4){
			cmap1_3d[ipx*(NY + 1)*(NZ + 1) + ipy*(NZ + 1) + ipz] = i;//test
		}
	}

	/*
	if (i < ncell){
		unsigned int state = cstate[i];
		cell_pos_set cme = cs[i];
		const bool afm_cell = state == ALIVE || state == FIX || state == MUSUME;
		int4 lat = {
			((int)(cme.x*inv_dx + NX)) % NX,
			((int)(cme.y*inv_dy + NY)) % NY,
			((int)(cme.z*inv_dz + NZ)) % NZ,
			0
		};


		if (!afm_cell){
			
			float rad = get_radius(state);
			for (int k = lat.x - irx; k <= lat.x + irx; k++){
				float mx = k*dx;
				int ipx = (k + NX) % NX;
				float diffx = p_diff_x(mx, cme.x);
				float diffxSq = diffx*diffx;

				for (int l = lat.y - iry; l <= lat.y + iry; l++){
					float my = l*dy;
					int ipy = (l + NY) % NY;
					float diffy = p_diff_y(my, cme.y);
					float diffySq = diffy*diffy;
					for (int ipz = zmin; ipz <= zmax; ipz++){
						float mz = ipz*dz;
						float diffz = mz - cme.z;
						if (diffz*diffz + diffySq + diffxSq < rad*rad){
							cmap1_3d[ipx*(NY + 1)*(NZ + 1) + ipy*(NZ + 1) + ipz] = i;
						}
					}
				}
			}
		}

		else{
			int _irz = 2 * irz;
			int _iry = 2 * iry;
			int _irx = 2 * irx;
			const int z_b = lat.z - _irz;
			const int z_u = lat.z + _irz;
			const int zmin = z_b >= 1 ? z_b : 1;
			const int zmax = z_u < (int)NZ ? z_u : (int)NZ - 1;
			float rad = get_radius(state);
			for (int k = lat.x - _irx; k <= lat.x + _irx; k++){
				float mx = k*dx;
				int ipx = (k + NX) % NX;
				float diffx = p_diff_x(mx, cme.x);
				float diffxSq = diffx*diffx;

				for (int l = lat.y - _iry; l <= lat.y + _iry; l++){
					float my = l*dy;
					int ipy = (l + NY) % NY;
					float diffy = p_diff_y(my, cme.y);
					float diffySq = diffy*diffy;
					for (int ipz = zmin; ipz <= zmax; ipz++){
						float mz = ipz*dz;
						float diffz = mz - cme.z;
						float distSq = diffz*diffz + diffySq + diffxSq;
						if (distSq < 4 * rad*rad){
							cmap2_3d[ipx*(NY + 1)*(NZ + 1) + ipy*(NZ + 1) + ipz] = 1.0f;
							if (distSq < rad*rad){
								cmap1_3d[ipx*(NY + 1)*(NZ + 1) + ipy*(NZ + 1) + ipz] = i;
							}
						}
					}
				}
			}
		}
		

	}
	*/
}

__global__ void setup_map_impl3(int ncell, cell_pos_set* cs, unsigned int* cstate, int* cmap1_3d, float* cmap2_3d){
	int i = blockIdx.x*blockDim.x + threadIdx.x;

	if (i < ncell){
		unsigned int state = cstate[i];
		cell_pos_set cme = cs[i];
		const bool afm_cell = state == ALIVE || state == FIX || state == MUSUME;
		int4 lat = {
			((int)(cme.x*inv_dx + NX)) % NX,
			((int)(cme.y*inv_dy + NY)) % NY,
			((int)(cme.z*inv_dz + NZ)) % NZ,
			0
		};

		int mirx = (int)(get_radius(state)*NX / LX);
		int miry = (int)(get_radius(state)*NY / LY);
		int mirz = (int)(get_radius(state)*NZ / LZ);

		if (!afm_cell){
			const int z_b = lat.z - mirz;
			const int z_u = lat.z + mirz;
			const int zmin = z_b >= 1 ? z_b : 1;
			const int zmax = z_u < (int)NZ ? z_u : (int)NZ - 1;
			float rad = get_radius(state);
			for (int k = lat.x - mirx; k <= lat.x + mirx; k++){
				int ipx = (k + NX) % NX;
				for (int l = lat.y - miry; l <= lat.y + miry; l++){
					int ipy = (l + NY) % NY;
					for (int ipz = zmin; ipz <= zmax; ipz++){
						cmap1_3d[ipx*(NY + 1)*(NZ + 1) + ipy*(NZ + 1) + ipz] = i;
					}
				}
			}
		}

		else{
			int _irz = 2 * mirz;
			int _iry = 2 * miry;
			int _irx = 2 * mirx;
			const int z_b = lat.z - _irz;
			const int z_u = lat.z + _irz;
			const int zmin = z_b >= 1 ? z_b : 1;
			const int zmax = z_u < (int)NZ ? z_u : (int)NZ - 1;
			float rad = get_radius(state);
			for (int k = lat.x - _irx; k <= lat.x + _irx; k++){
				int ipx = (k + NX) % NX;
				int ipx2 = (lat.x + (k - lat.x) / 2 + NX) % NX;
				for (int l = lat.y - _iry; l <= lat.y + _iry; l++){
					int ipy = (l + NY) % NY;
					int ipy2 = (lat.y + (l - lat.y) / 2 + NY) % NY;
					for (int ipz = zmin; ipz <= zmax; ipz++){
						int ipz2 = lat.z + (ipz - lat.z) / 2;
							cmap2_3d[ipx*(NY + 1)*(NZ + 1) + ipy*(NZ + 1) + ipz] = 1.0f;

							cmap1_3d[ipx2*(NY + 1)*(NZ + 1) + ipy2*(NZ + 1) + ipz2] = i;
						
					}
				}
			}
		}

	}
}

__global__ void setup_map_impl4(int ncell, cell_pos_set* cs, unsigned int* cstate, int* cmap1_3d, float* cmap2_3d){
	int i = blockIdx.x*blockDim.x + threadIdx.x;

	if (i < ncell){
		__shared__ char4 s_map_memb_sphere[MEMB_SPHERE_COUNT]; __shared__ int s_map_memb_count;
		__shared__ char4 s_map_normal_sphere[NSPHERE_COUNT]; __shared__ int s_map_normal_count;
		__shared__ char4 s_map_twice_sphere[T_SPHERE_COUNT]; __shared__ int s_map_twice_count;
		if (threadIdx.x == 0){
			s_map_memb_count = map_memb_count;
			s_map_normal_count = map_normal_count;
			s_map_twice_count = map_twice_count;
			memcpy(s_map_memb_sphere, map_memb_sphere, sizeof(char4) * MEMB_SPHERE_COUNT);
			memcpy(s_map_normal_sphere, map_normal_sphere, sizeof(char4) * NSPHERE_COUNT);
			memcpy(s_map_twice_sphere, map_twice_sphere, sizeof(char4) * T_SPHERE_COUNT);
		}
		unsigned int state = cstate[i];
		cell_pos_set cme = cs[i];
		const bool afm_cell = state == ALIVE || state == FIX || state == MUSUME;
		int4 lat = {
			(int)round(cme.x*inv_dx),
			(int)round(cme.y*inv_dy),
			(int)round(cme.z*inv_dz),
			0
		};
		__syncthreads();
		
		if (state == MEMB){
			for (int j = 0; j < s_map_memb_count; j++){
				int ipx = (lat.x + s_map_memb_sphere[j].x + NX) % NX;
				int ipy = (lat.y + s_map_memb_sphere[j].y + NY) % NY;
				int ipz = lat.z + s_map_memb_sphere[j].z;

				if (ipz >= 0 && ipz < NZ)cmap1_3d[ipx*(NY + 1)*(NZ + 1) + ipy*(NZ + 1) + ipz] = i;
			}
		}
		else{
			for (int j = 0; j < s_map_normal_count; j++){
				int ipx = (lat.x + s_map_normal_sphere[j].x + NX) % NX;
				int ipy = (lat.y + s_map_normal_sphere[j].y + NY) % NY;
				int ipz = lat.z + s_map_normal_sphere[j].z;

				if (ipz >= 0 && ipz < NZ)cmap1_3d[ipx*(NY + 1)*(NZ + 1) + ipy*(NZ + 1) + ipz] = i;
			}
		}
		if (afm_cell){

			for (int j = 0; j < s_map_twice_count; j++){
				int ipx = (lat.x + s_map_twice_sphere[j].x + NX) % NX;
				int ipy = (lat.y + s_map_twice_sphere[j].y + NY) % NY;
				int ipz = lat.z + s_map_twice_sphere[j].z;

				if (ipz >= 0 && ipz < NZ)cmap2_3d[ipx*(NY + 1)*(NZ + 1) + ipy*(NZ + 1) + ipz] = 1.0f;
			}
		}
		

	}
}
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
	cudaMemcpy(map_normal_sphere, nm, sizeof(char4) * NSPHERE_COUNT,cudaMemcpyHostToDevice);
	cudaMemcpy(&map_normal_count, &count, sizeof(int), cudaMemcpyHostToDevice);

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

	cudaMemcpy(map_memb_sphere, nm, sizeof(char4) * MEMB_SPHERE_COUNT, cudaMemcpyHostToDevice);
	cudaMemcpy(&map_memb_count, &count, sizeof(int), cudaMemcpyHostToDevice);

	count = 0;
	for (char i = -t_irx; i <= t_irx; i++){
		float x = dx*i;
		for (char j = -t_iry; j <= t_iry; j++){
			float y = dy*j;
			for (char k = -t_irz; k <= t_irz; k++){
				float z = dz*k;
				if (x*x + y*y + z*z < 4*R_max*R_max){
					nm[count].x = i;
					nm[count].y = j;
					nm[count].z = k;
					count++;
				}
			}
		}
	}

	cudaMemcpy(map_twice_sphere, nm, sizeof(char4) * T_SPHERE_COUNT, cudaMemcpyHostToDevice);
	cudaMemcpy(&map_twice_count, &count, sizeof(int), cudaMemcpyHostToDevice);

}
void setup_map(DeviceData*d, int* cmap1_3d, float* cmap2_3d){
	cudaMemset(cmap1_3d, 0x80, sizeof(int)*(NX + 1)*(NY + 1)*(NZ + 1));
	cudaMemset(cmap2_3d, 0, sizeof(float)*(NX + 1)*(NY + 1)*(NZ + 1));

	setup_map_impl4<<<d->ncell/256+1,256>>>(d->ncell, d->c_pos_d[d->current], (unsigned int*)d->c_state_d, cmap1_3d, cmap2_3d);
}