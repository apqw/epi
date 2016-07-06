#include "global_cuda.h"
#include "ext_stim.h"
__constant__ static const float kb = 0.025f;
__constant__ static const float DUR_ALIVE = 0.5f;
__constant__ static const float DUR_DEAD = 2.0f;
__constant__ static const float DB = 0.0009f;
__device__  float fB(float age, float B, bool cornif) {


	using namespace cont;
	return (cornif&&age > THRESH_DEAD - DUR_ALIVE&&age <= THRESH_DEAD + DUR_DEAD ? 1.0f : 0.0f) - kb*B;
}
#define MIDX(x,y,z) ((x)*(NY+1)*(NZ+1)+(y)*(NZ+1)+(z))
__global__ void calc_ext_stim_impl(float* c_ext_stim, float* new_ext_stim, int* cmap1, float* cmap2, unsigned int* cstate, float* agek_arr){
	
	//const int iz_bound = (int)((zzmax + 2.0*R_max) *inv_dz);
	const int ix = blockIdx.x;
	const int iy = blockIdx.y;
	const int iz = threadIdx.x;

	const int next_ix = (ix + 1) % NX;
	const int prev_ix = (ix - 1 + NX) % NX;

	const int next_iy = (iy + 1) % NY;
	const int prev_iy = (iy - 1 + NY) % NY;

	const int next_iz = iz == NZ - 1 ? NZ - 2 : iz + 1;
	const int prev_iz = iz == 0 ? 1 : iz - 1;

	const int center = MIDX(ix, iy, iz);
	float dum_age = 0.0f;
	bool flg_cornified = false;
	
	if (cmap1[center] >= 0){
		unsigned int state = cstate[cmap1[center]];
		
		if (state == DEAD || state == ALIVE){
			dum_age = agek_arr[cmap1[center]];
			flg_cornified = true;
		}
	}
	
	/*
	auto& cext = carr[j][k][l];

	narr[j][k][l] = cext
		+ DT_Ca*(DB *
		(cmap2()[prev_x][k][l] * (carr[prev_x][k][l] - cext)
		+ cmap2()[j][prev_y][l] * (carr[j][prev_y][l] - cext)
		+ cmap2()[j][k][prev_z] * (carr[j][k][prev_z] - cext)
		+ cmap2()[j][k][next_z] * (carr[j][k][next_z] - cext)
		+ cmap2()[j][next_y][l] * (carr[j][next_y][l] - cext)
		+ cmap2()[next_x][k][l] * (carr[next_x][k][l] - cext)

		) *inv_dz*inv_dz
		+ fB(dum_age, cext, flg_cornified));
		*/
	using namespace cont;
	float cvalue = c_ext_stim[center];
	//float a = new_ext_stim[center];
	
	new_ext_stim[center] = cvalue + DT_Ca*(DB*(
		cmap2[MIDX(prev_ix, iy, iz)] * (c_ext_stim[MIDX(prev_ix, iy, iz)] - cvalue)
		+ cmap2[MIDX(ix, prev_iy, iz)] * (c_ext_stim[MIDX(ix, prev_iy, iz)] - cvalue)
		+ cmap2[MIDX(ix, iy, prev_iz)] * (c_ext_stim[MIDX(ix, iy, prev_iz)] - cvalue)
		+ cmap2[MIDX(ix, iy, next_iz)] * (c_ext_stim[MIDX(ix, iy, next_iz)] - cvalue)
		+ cmap2[MIDX(ix, next_iy, iz)] * (c_ext_stim[MIDX(ix, next_iy, iz)] - cvalue)
		
	
		
		+ cmap2[MIDX(next_ix, iy, iz)] * (c_ext_stim[MIDX(next_ix, iy, iz)] - cvalue)
		
			
			
		)*inv_dz*inv_dz + fB(dum_age, cvalue, flg_cornified));
			return;
			//);
	
	
}

void calc_ext_stim(DeviceData*d,float* c_ext_stim, float* new_ext_stim, int* cmap1, float* cmap2){
	calc_ext_stim_impl << < dim3(NX, NY), NZ >> >(c_ext_stim, new_ext_stim, cmap1, cmap2, (unsigned int*)d->c_state_d, d->c_agek_d);
}