#include "define.h"
#include "CData.h"
#include "CellDataMan.h"
#define PRE_GRID_SIZE (R_max*CDEF(2.0))

#define CALC_GRID_DIV(x) ((int)((x)/PRE_GRID_SIZE)-1)
#define ANX CALC_GRID_DIV(LX)
#define ANY CALC_GRID_DIV(LY)
#define ANZ CALC_GRID_DIV(LZ)

#define CALC_GRID_POS_X(x) ((int)(ANX*(x)/LX))
#define CALC_GRID_POS_Y(y) ((int)(ANY*(y)/LY))
#define CALC_GRID_POS_Z(z) ((int)(ANZ*(z)/LZ))



#define MEMB_FOLD_LAYER_ASSUMPTION (2)

#define GRID_STORE_MAX (2*2*2+(int)(2*COMPRESS_FACTOR*2*COMPRESS_FACTOR*MEMB_FOLD_LAYER_ASSUMPTION))

struct CellIndexStore{
	CellIndex arr[GRID_STORE_MAX];
};

using CellIndexGrid = CArr3D<CellIndexStore, ANX, ANY, ANZ>;
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
__global__ void build_grid(const CValue<int>::ptr_type ncell, const CPosArr::ptr_type current_pos, CellIndexGrid::ptr_type idx_grd, CellCountGrid::ptr_type count_grd){
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
		const int my_count = atomicAdd(count_grd.get() + gpos, 1);
		thrust::raw_reference_cast(idx_grd[gpos]).arr[my_count] = index;
		assert(my_count < GRID_STORE_MAX);
	}
}

void reset_grid_count(CellCountGrid* cg){
	cg->_memset(0);
}

__global__ void connect_proc(const int ncell, const CellIndex start
	, const CPosArr::ptr_type cpos, const CArr<CELL_STATE>::ptr_type cstate
	, CArr<CellConnectionData>::ptr_type conndata
	, const CellIndexGrid::ptr_type idx_grd, const CellCountGrid::ptr_type count_grd){
	
	const int my_proc_id = threadIdx.x % 32;
	if (my_proc_id >= 3 * 3 * 3)return;
	
	
	//__shared__ CellIndex start;
	__shared__ CellPos pos_set[32];
	__shared__ int4 grd_set[32];

	const int virtual_id = threadIdx.x / 32;
	__syncthreads();
	
	
	
	const CellIndex index = blockDim.x*blockIdx.x / 32 + virtual_id + start;
	if (index < ncell){
		
		if (my_proc_id == 0){

			pos_set[virtual_id] = cpos[index];
			grd_set[virtual_id] = make_int4(
				CALC_GRID_POS_X(pos_set[virtual_id].x),
				CALC_GRID_POS_Y(pos_set[virtual_id].y),
				CALC_GRID_POS_Z(pos_set[virtual_id].z), 0);
		}
		__syncthreads();
		const CellPos my_pos = pos_set[virtual_id];
		const int4 my_grd = grd_set[virtual_id];
		
		if (my_proc_id == 0){
#ifdef DBG
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
#endif
		}
		const int proc_x = my_proc_id%3+my_grd.x-1;
		const int proc_y = (my_proc_id / 3) % 3 + my_grd.y - 1;
		const int proc_z = (my_proc_id / 9) +my_grd.z - 1;
		const int arr_idx = gi::idx(proc_x, proc_y, proc_z);
		const int stored_num = count_grd[arr_idx];
		const CellIndexStore& conn = thrust::raw_reference_cast(idx_grd[arr_idx]);
		for (int i = 0; i < stored_num; i++){
			const CellIndex opidx = conn.arr[i];
			if (index <= opidx)continue;
			const CellPos oppos = cpos[opidx];
			const CELL_STATE cst = cstate[opidx];
			const real rad_sum = R_max + cst == MEMB ? R_memb : R_max;
			if (p_dist_sq(my_pos, oppos) <= LJ_THRESH*LJ_THRESH*rad_sum*rad_sum){
				CellConnectionData& my_cdata = thrust::raw_reference_cast(conndata[index]);
				CellConnectionData& op_cdata = thrust::raw_reference_cast(conndata[opidx]);
				my_cdata.connect_index[atomicAdd(&my_cdata.connect_num, 1)] = opidx;
				op_cdata.connect_index[atomicAdd(&op_cdata.connect_num, 1)] = index;
			}
		}


	}

}

void connect_cell(CellDataMan* cm){
	static CellIndexGrid cidxg;
	static CellCountGrid ccntg;

	reset_grid_count(&ccntg);

	reset_connect_count << <DEFAULT_THB_ALL_CELL >> >(cm->ncell, cm->nmemb, cm->connection_data);
	build_grid << <DEFAULT_THB_ALL_CELL >> >(cm->ncell, cm->pos.current(), cidxg, ccntg);
	gpuErrchk(cudaDeviceSynchronize());

	const CellIndex start = (int)cm->nmemb + (int)cm->nder;
	const int ncell = cm->ncell;
	connect_proc << <(ncell - start - 1) / 32 + 1, 32 * 32 >> >(ncell, start, cm->pos.current(), cm->state, cm->connection_data, cidxg, ccntg);
	gpuErrchk(cudaDeviceSynchronize());
	//need connect reset


}