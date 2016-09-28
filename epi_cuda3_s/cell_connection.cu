#include "define.h"
#include "CData.h"
#include "CellDataMan.h"
#include <cuda_fp16.h>


#define CELL_IN_BLOCK 2
#define PAR_NUM 8
struct half3{
	half x;
	half y;
	half z;
};
/*
struct CellPackedInfo{
	half x;
	half y;
	half z;
	unsigned short count;

	CellIndex idx;
	__device__ void set_info(CellPos cp,CellIndex i){
		x = __float2half(cp.x);
		y = __float2half(cp.y);
		z = __float2half(cp.z);
		idx = i;
	}
	__device__ float3 get_pos()const{
		return make_float3(__half2float(x), __half2float(y), __half2float(z));
	}
};
struct CellInfoStore{
	int count;
	CellPackedInfo arr[GRID_STORE_MAX];
	__device__ int get_stored_num()const{
		return count;
	}

	__device__ int* get_stored_num_ptr(){
		return &count;
	}
};

using CellInfoGrid = CArr3D<CellInfoStore, ANX, ANY, ANZ>;
*/
using CellCountGrid = CArr3D<int, ANX, ANY, ANZ>;
using gi = __midx<ANX, ANY, ANZ>;

void grid_condition_check(){
	static_assert(ANX >= 3, "X grid div num has to be >=3");
	static_assert(ANY >= 3, "Y grid div num has to be >=3");
	static_assert(ANZ >= 3, "Z grid div num has to be >=3");
}

__global__ void reset_connect_count(CValue<int>::ptr_type ncell,CValue<int>::ptr_type nmemb, CArr<CellConnectionData>::ptr_type conn){
	const CellIndex index = blockDim.x*blockIdx.x + threadIdx.x;
	if (index < (int)ncell[0]){
		thrust::raw_reference_cast(conn[index]).connect_num = index < (int)nmemb[0] ? MEMB_CONN_NUM : 0;
	}
}
__global__ void build_grid(const CValue<int>::ptr_type ncell,const CValue<int>::ptr_type nmemb, const CPosArr::ptr_type current_pos, CellInfoGrid::ptr_type idx_grd){
	const CellIndex index = blockDim.x*blockIdx.x + threadIdx.x;
	if (index < (int)ncell[0]){
		const CellPos cpos = current_pos[index];
		const int aix = CALC_GRID_POS_X(cpos.x);
		const int aiy = CALC_GRID_POS_Y(cpos.y);
		const int aiz = CALC_GRID_POS_Z(cpos.z);

		if (!(aix < ANX&&aix >= 0)){
			printf("aix error:%d\n", aix);
			assert(aix < ANX&&aix >= 0);
		}
		if (!(aiy < ANY&&aiy >= 0)){
			printf("aiy error:%d\n", aiy);
			assert(aiy < ANY&&aiy >= 0);
		}
		if (!(aiz < ANZ&&aiz >= 0)){
			printf("aiz error:%d\n", aiz);
			assert(aiz < ANZ&&aiz >= 0);
		}
		const int gpos = gi::idx(aix, aiy, aiz);
		const int my_count = atomicAdd(thrust::raw_reference_cast(idx_grd[gpos]).get_stored_num_ptr(), 1);
		/*
		printf("%d count\n", my_count);
		if (!(my_count < GRID_STORE_MAX)){
			for (int i = 0; i < 100; i++){
				printf("broken conn count:%d\n", my_count);
			}

			assert(my_count < GRID_STORE_MAX);
		}
		*/
		assert(my_count < GRID_STORE_MAX);
		thrust::raw_reference_cast(idx_grd[gpos]).arr[my_count] = index;
		

		//if (thrust::raw_reference_cast(idx_grd[gpos]).arr[my_count].idx >= MAX_CELL_NUM)printf("%d broke idx\n", thrust::raw_reference_cast(idx_grd[gpos]).arr[my_count].idx);
	}
}

void reset_grid_count(CellInfoGrid* cg){
	cg->_memset(0);
}

__global__ void connect_proc(const int ncell,const int nmemb
	, const CPosArr cpos
	, CArr<CellConnectionData>::ptr_type conndata
	, const CellInfoGrid::ptr_type idx_grd){
	
	const int my_proc_id = (threadIdx.x/PAR_NUM) % (32);
	const int my_par_offset = threadIdx.x%PAR_NUM;
	if (my_proc_id >= 3 * 3 * 3)return;
	
	
	//__shared__ CellIndex start;
	__shared__ real3 pos_set[CELL_IN_BLOCK];
	__shared__ int3 grd_set[CELL_IN_BLOCK];
	__shared__ CellConnectionData* myconn[CELL_IN_BLOCK];

	const int virtual_id = threadIdx.x / (PAR_NUM * 32);
	__syncthreads();
	
	
	
	const CellIndex index = blockIdx.x*CELL_IN_BLOCK + virtual_id + nmemb;
	if (index < ncell){
		
		if (my_proc_id == 0&&my_par_offset==0){
			const CellPos tmpc = cpos[index];
			pos_set[virtual_id] = make_real3(tmpc.x,tmpc.y,tmpc.z);
			grd_set[virtual_id] = make_int3(
				CALC_GRID_POS_X(pos_set[virtual_id].x),
				CALC_GRID_POS_Y(pos_set[virtual_id].y),
				CALC_GRID_POS_Z(pos_set[virtual_id].z));
			myconn[virtual_id]=thrust::raw_pointer_cast(conndata + index);
		}
		__syncthreads();
		const real3& my_pos = pos_set[virtual_id];
		const int3& my_grd = grd_set[virtual_id];
#ifdef DBG
		if (my_proc_id == 0){

			if (!(my_grd.x < ANX&&my_grd.x >= 0)){
				printf("my aix error:%d\n", my_grd.x);
				assert((my_grd.x < ANX&&my_grd.x >= 0));
			}
			if (!(my_grd.y < ANY&&my_grd.y >= 0)){
				printf("my aiy error:%d\n", my_grd.y);
				assert((my_grd.y < ANY&&my_grd.y >= 0));
			}
			if (!(my_grd.z < ANZ&&my_grd.z >= 0)){
				printf("my aiz error:%d\n", my_grd.z);
				assert((my_grd.z < ANZ&&my_grd.z >= 0));
			}
			

		}
#endif
		const int proc_x = my_proc_id%3+my_grd.x-1;
		const int proc_y = (my_proc_id / 3) % 3 + my_grd.y - 1;
		const int proc_z = (my_proc_id / 9) +my_grd.z - 1;
		const int arr_idx = gi::idx(proc_x, proc_y, proc_z);
		printf("%d arridx\n", arr_idx);
		const CellInfoStore& conn = thrust::raw_reference_cast(idx_grd[arr_idx]);

		const int stored_num = conn.count;
		//if(threadIdx.x%10==0)printf("%d aa\n", stored_num);
		for (int i = my_par_offset; i < stored_num; i+=PAR_NUM){
			
			const CellIndex opidx = conn.arr[i];
			if (index <= opidx)continue;
			const CellPos oppos = cpos[opidx];
			const real rad_sum = R_max + (opidx<nmemb ? R_memb : R_max);
			const real diffx = p_diff_x(my_pos.x, oppos.x);
			const real diffy = p_diff_y(my_pos.y, oppos.y);
			const real diffz = my_pos.z- oppos.z;
			
			if (diffx*diffx+diffy*diffy+diffz*diffz <= LJ_THRESH*LJ_THRESH*rad_sum*rad_sum){
				
				
				CellConnectionData* op_cdata = thrust::raw_pointer_cast(conndata + opidx);
				
				myconn[virtual_id]->connect_index[atomicAdd(&myconn[virtual_id]->connect_num, 1)] = opidx;
				
				//if(opidx>=MAX_CELL_NUM)printf("%d opidx\n", opidx);
				op_cdata->connect_index[atomicAdd(&op_cdata->connect_num, 1)] = index;
				
			}
			
		}


	}

}

__global__ void refresh_gj(const CValue<int>::ptr_type ncell,CArr<CellConnectionData>::ptr_type conn){
	const CellIndex index = blockDim.x*blockIdx.x + threadIdx.x;
	if (index < (int)ncell[0]){
		auto& connref = thrust::raw_reference_cast(conn[index]);
		using HashType = decltype(connref.gj);
		connref.gj.foreach([]__host__ __device__(HashType::node_type*node, CellIndex key, gj_value& value){
			value.checked = false;
		});

		const int cnum = connref.connect_num;
		for (int i = 0; i < cnum; i++){
			connref.gj.at_with_emplace(connref.connect_index[i]).checked = true;
		}

		connref.gj.foreach([]__host__ __device__(HashType::node_type*node, CellIndex key, gj_value& value){
			if (!value.checked){
				value.value = gj_init;
			}
		});
	}
}

__global__ void see_connection(CArr<CellConnectionData>::ptr_type conn){
	const CellIndex index = blockDim.x*blockIdx.x + threadIdx.x;
	if (index < 60000){
		CellConnectionData& my_cdata = thrust::raw_reference_cast(conn[index]);

		if(index%100==0)printf("conn num:%d\n", my_cdata.connect_num);
	}
}

void connect_cell(CellDataMan* cm){
	//static CellInfoGrid cidxg;
	//static CellCountGrid ccntg;
	//for test
	reset_grid_count(&cm->cell_info_grid);
	reset_connect_count << <DEFAULT_THB_ALL_CELL >> >(cm->ncell, cm->nmemb, cm->connection_data);
	build_grid << <DEFAULT_THB_ALL_CELL >> >(cm->ncell, cm->nmemb, cm->pos.current(), cm->cell_info_grid);
	gpuErrchk(cudaDeviceSynchronize());
	const int nmemb = cm->nmemb;
	const int ncell = cm->ncell;
	connect_proc << <(ncell - nmemb - 1) / CELL_IN_BLOCK + 1, PAR_NUM*CELL_IN_BLOCK * 32 >> >(ncell, nmemb, cm->pos.current(), cm->connection_data, cm->cell_info_grid);
	gpuErrchk(cudaDeviceSynchronize());
	refresh_gj << <DEFAULT_THB_ALL_CELL >> >(cm->ncell, cm->connection_data);
	gpuErrchk(cudaDeviceSynchronize());
	//cudaUnbindTexture(pos_tex);
	/*
	see_connection << <DEFAULT_THB_ALL_CELL >> >(cm->connection_data);
	gpuErrchk(cudaDeviceSynchronize());
	*/
	//need connect reset


}
