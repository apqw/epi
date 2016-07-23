#pragma once
#include <thrust/device_ptr.h>
#include <string>
#include "define.h"
#include "CData.h"
#include "switcher.h"
#include "hashmap.h"
#include <cuda_fp16.h>
using RField3DSw = Field3DMulti<real, 2>;
using CPosArrSw = CArrMulti<CellPos, MAX_CELL_NUM, 2>;
using CPosArr = CPosArrSw::Arr_type;
using RArrSw = CArrMulti<real, MAX_CELL_NUM, 2>;

using FieldMask3D = Field3D<FieldMask_t>;
using IndexMap3D = Field3D<CellIndex>;

struct CellConnectionData{
	int connect_num;
	CellIndex connect_index[MAX_CONNECT_CELL_NUM];
	IntegerHashmap<real> gj;
	
	//CellConnectionData(){}
};

struct CellPackedInfo{
	half x;
	half y;
	half z;
	unsigned short count;

	CellIndex idx;
	__device__ void set_info(CellPos cp, CellIndex i);
	__device__ float3 get_pos()const;
};
struct CellInfoStore{
	int count;
	CellIndex arr[GRID_STORE_MAX];
	__device__ int get_stored_num()const{
		return count;
	}

	__device__ int* get_stored_num_ptr(){
		return &count;
	}
};
using CellInfoGrid = CArr3D<CellInfoStore, ANX, ANY, ANZ>;

//__global__ void cell_connection_device_init(CArr<CellConnectionData> aa);
void reset_index_map(IndexMap3D* imap);

void reset_field_mask(FieldMask3D* fmask);
class CellDataMan
{
	CPosArrSw _pos;
	RArrSw _ca2p;
	RArrSw _IP3;
	RField3DSw _ATP;
	RField3DSw _ext_stim;

	switcher swt;
	switcher ca2p_swt;
	void _load_from_file(const std::string& path);
	__device__ __host__ void _set_cell(
		CellIndex index,
		CELL_STATE state,
		int fix_origin,
		real x,real y,real z,
		//rad,
		real ca2p,
		real ca2p_avg,
		real IP3,
		real ex_inert,
		real agek,
		real ageb,
		real ex_fat,
		real fat,
		real spr_len,
		int div_times,
		bool is_malig
		);
public:
	CArr<CELL_STATE> state;

	CArr<int> fix_origin;
	CArr<CellConnectionData> connection_data;
	CArr<CellIndex> pair_index;
	
	CArr<real> ca2p_avg;
	CArr<real> ex_inert;

	CArr<real> agek;
	CArr<real> ageb;
	CArr<real> ex_fat;
	CArr<real> in_fat;
	CArr<real> spr_nat_len;
	//CArr<int> uid;
	CArr<int> rest_div_times;
	CArr<CellIndex> dermis_index;
	CArr<bool> is_malignant;

	CValue<real> zmax; //one value
	CValue<int> ncell;
	CValue<int> nder;
	CValue<int> nmemb;
	CValue<int> sw;
	//CValue<int> next_uid;

	SwtAccessor<decltype(_ATP)> ATP;
	SwtAccessor<decltype(_ca2p)> ca2p;
	SwtAccessor<decltype(_IP3)> IP3;
	SwtAccessor<decltype(_pos)> pos;
	SwtAccessor<decltype(_ext_stim)> ext_stim;

	FieldMask3D field_mask;
	IndexMap3D index_map;
	CellInfoGrid cell_info_grid;
	
	CellDataMan();
	~CellDataMan();

	__host__ void cell_phase_switch();
	__host__ void ca2p_phase_switch();

	bool check_initialized()const;
	void init(const std::string& init_data_path, bool use_last = false, const std::string& init_uvp_data = "", const std::string& init_w_data = "", const std::string& init_ATP_data="", const std::string& init_ext_stim_data="");
	
};

