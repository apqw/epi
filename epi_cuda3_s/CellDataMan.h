#pragma once
#include <thrust/device_ptr.h>
#include <string>
#include "define.h"
#include "CData.h"
#include "switcher.h"

using RField3DSw = Field3DMulti<real, 2>;
using CPosArrSw = CArrMulti<CellPos, MAX_CELL_NUM, 2>;
using RArrSw = CArrMulti<real, MAX_CELL_NUM, 2>;

using FieldMask3D = Field3D<float>;
using IndexMap3D = Field3D<CellIndex>;
struct CellConnectionData{
	CArr<CellIndex> connect_index;
	CArr<real> gj;
	CellConnectionData();
};


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
	CArr<int> uid;
	CArr<int> rest_div_times;
	CArr<CellIndex> dermis_index;
	CArr<bool> is_malignant;

	CValue<real> zmax; //one value
	CValue<int> ncell;
	CValue<int> nder;
	CValue<int> nmemb;
	CValue<int> sw;
	CValue<int> next_uid;

	SwtAccessor<decltype(_ATP)> ATP;
	SwtAccessor<decltype(_ca2p)> ca2p;
	SwtAccessor<decltype(_IP3)> IP3;
	SwtAccessor<decltype(_pos)> pos;
	SwtAccessor<decltype(_ext_stim)> ext_stim;

	FieldMask3D field_mask;
	IndexMap3D index_map;
	
	CellDataMan();
	~CellDataMan();

	__host__ void cell_phase_switch();
	__host__ void ca2p_phase_switch();

	bool check_initialized()const;
	void init(const std::string& init_data_path, bool use_last = false, const std::string& init_uvp_data = "", const std::string& init_w_data = "");
	
};

