
#include "define.h"
#include "utils.h"
#include "CellManager.h"
#define kb (0.025f)
#define DUR_ALIVE (0.5f)
#define DUR_DEAD (2.0f)
#define DB (0.0009f)
__device__  real fB(real age, real B, bool cornif) {



	return (cornif&&age > THRESH_DEAD - DUR_ALIVE&&age <= THRESH_DEAD + DUR_DEAD ? 1.0f : 0.0f) - kb*B;
}

/*
	cache 3 planes
*/
__global__ void calc_ext_stim_impl(const real* __restrict__ c_ext_stim, real* __restrict__ new_ext_stim, const int* __restrict__ cmap1, const float* __restrict__ cmap2,
	const unsigned int* __restrict__ cstate, const real* __restrict__ agek_arr, int iz_bound, int shmem_border){

	//const int iz_bound = (int)((zzmax + 2.0*R_max) *inv_dz);
	const int iz = threadIdx.x;
	const int ix = blockIdx.x;
	const int iy = blockIdx.y;
	//if (iz > NZ||iy>NY||ix>NX)return;
	extern __shared__ float cmap2_line[];
	real* ext_stim_line = (real*)&cmap2_line[shmem_border];
	//__shared__ real ext_stim_line[NZ + 1];
	const int next_iy = (iy + 1) % NY;
	const int prev_iy = (iy - 1 + NY) % NY;


	const int next_ix = (ix + 1) % NX;
	const int prev_ix = (ix - 1 + NX) % NX;

	const int next_iz = iz == NZ - 1 ? NZ - 2 : iz + 1;
	const int prev_iz = iz == 0 ? 1 : iz - 1;

	const int center = midx<NX + 1, NY + 1, NZ + 1>(ix, iy, iz);
	
	const int idx1 = midx<NX + 1, NY + 1, NZ + 1>(ix, prev_iy, iz);
	const int idx2 = midx<NX + 1, NY + 1, NZ + 1>(ix, next_iy, iz);
	const int idx3 = midx<NX + 1, NY + 1, NZ + 1>(prev_ix, iy, iz);
	const int idx4 = midx<NX + 1, NY + 1, NZ + 1>(next_ix, iy, iz);
	cmap2_line[iz] = cmap2[center];
	const real cvalue= (ext_stim_line[iz] = c_ext_stim[center]);
	if (iz == iz_bound - 1){
		cmap2_line[next_iz] = cmap2[midx<NX + 1, NY + 1, NZ + 1>(ix, iy, next_iz)];
		ext_stim_line[next_iz] = c_ext_stim[midx<NX + 1, NY + 1, NZ + 1>(ix, iy, next_iz)];
	}

	if (iz == 0){
		cmap2_line[prev_iz] = cmap2[midx<NX + 1, NY + 1, NZ + 1>(ix, iy, prev_iz)];
		ext_stim_line[prev_iz] = c_ext_stim[midx<NX + 1, NY + 1, NZ + 1>(ix, iy, prev_iz)];
	}
	




	
	
	real dum_age = 0.0f;
	bool flg_cornified = false;

	if (cmap1[center] >= 0){
		unsigned int state = cstate[cmap1[center]];

		if (state == DEAD || state == ALIVE){
			dum_age = agek_arr[cmap1[center]];
			flg_cornified = true;
		}
	}
	__syncthreads();
	new_ext_stim[center] = cvalue + DT_Ca*(DB*(
		cmap2_line[prev_iz] * (ext_stim_line[prev_iz] - cvalue)
		+ cmap2_line[next_iz] * (ext_stim_line[next_iz] - cvalue)
		+ cmap2[idx1] * (c_ext_stim[idx1] - cvalue)
		+ cmap2[idx2] * (c_ext_stim[idx2] - cvalue)
		+ cmap2[idx3] * (c_ext_stim[idx3] - cvalue)
		+ cmap2[idx4] * (c_ext_stim[idx4] - cvalue)
		)*inv_dz*inv_dz + fB(dum_age, cvalue, flg_cornified));
	
	return;
	//);


}

void calc_ext_stim(CellManager*cm, real* c_ext_stim, real* new_ext_stim, int* cmap1, float* cmap2, real zmax){
	int iz_bound=(int)((zmax + 2.0f*R_max) *inv_dz)+1;
	int shmem_border = iz_bound +1+ (iz_bound + 1) % 2;
	//reduce thread num -> faster
	calc_ext_stim_impl << < dim3(NX, NY), iz_bound, (shmem_border)*sizeof(float) + (iz_bound + 1)*sizeof(real) >> >(c_ext_stim, new_ext_stim, cmap1, cmap2, cm->state, cm->agek, iz_bound, shmem_border);
	cudaDeviceSynchronize();
}