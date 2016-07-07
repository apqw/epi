#include "global_cuda.h"
//#include <thrust\extrema.h>
//#include <thrust\device_ptr.h>
#include "cuda_host_util.h"
#define ZMAX_APPROX_PREC (1000.0f)
/*
struct opp{
	__host__ __device__ bool operator()(const cell_pos_set c1, const cell_pos_set c2){
		return c1.z < c2.z;
	};
};
float calc_zzmax(cell_pos_set* pos, int ncell){
	thrust::device_ptr<cell_pos_set> pos_dptr(pos+10);
	cell_pos_set temp;
	//*thrust::max_element(pos_dptr, pos_dptr + 1, opp())
	
	return pos_dptr[0].operator cell_pos_set().z;
}

__global__ void calc_zzmax_fool_impl(float* out, cell_pos_set* cs, int ncell){
	//extern __shared__ 
	int idx = 0;
	for (int i = 1; i < ncell; i++){
		if (cs[idx].z < cs[i].z)idx = i;
	}
	*out = cs[idx].z;
}

float calc_zzmax_2(DeviceData*d){
	calc_zzmax_fool_impl<<<1,1>>>(d->zzmax, d->c_pos_d[d->current], d->ncell);
	float mine;
	cudaMemcpy(&mine, d->zzmax, sizeof(float), cudaMemcpyDeviceToHost);
	return mine;
}
*/
__global__ void calc_zzmax_approx_pass1(int* out, cell_pos_set* cs, int ncell){
	__shared__ int thread_max;
	if (threadIdx.x == 0){
		thread_max = 0;
	}
	__syncthreads();
	int index = blockDim.x*blockIdx.x + threadIdx.x;
	if (index < ncell){
		int tmp = (int)(cs[index].z*ZMAX_APPROX_PREC);
		atomicMax(&thread_max, tmp);
		__syncthreads();
		if (threadIdx.x == 0){
			out[blockIdx.x] = thread_max;
		}
	}
}

__global__ void calc_zzmax_approx_pass2(float* out, int* block_max){
	__shared__ int final_max;
	
	if (threadIdx.x == 0){
		final_max = 0;
	}
	__syncthreads();
	atomicMax(&final_max, block_max[threadIdx.x]);
	__syncthreads();
	if (threadIdx.x == 0){
		*out = ((float)final_max) / ZMAX_APPROX_PREC;
	}
}

float calc_zzmax_approx(DeviceData*d){
	cudaMemset(d->block_max_store, 0, sizeof(int) * 512);
	//cudaMemset(d->zzmax, 0, sizeof(float));
	calc_zzmax_approx_pass1 << <512, 128 >> >(d->block_max_store, d->c_pos_d[d->current], d->ncell);
	cudaThreadSynchronize();
	calc_zzmax_approx_pass2 << <1, 512 >> >(d->zzmax, d->block_max_store);
	//cudaThreadSynchronize();
	float mine;
	cudaMemcpy(&mine, d->zzmax, sizeof(float), cudaMemcpyDeviceToHost); //sync
	return mine;
}