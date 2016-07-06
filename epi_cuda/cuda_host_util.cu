#include "global_cuda.h"
#include <thrust\extrema.h>
#include <thrust\device_ptr.h>
#include "cuda_host_util.h"
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