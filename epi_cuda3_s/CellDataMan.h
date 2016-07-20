#pragma once
#include <thrust/device_ptr.h>
#include "define.h"
#include "CData.h"
#include "switcher.h"

struct CellConnectionData{
	CArr<CellIndex> connect_index;
	CArr<real> gj;
	CellConnectionData();
};


class CellDataMan
{
	CArrMulti<CellPos, MAX_CELL_NUM,2> pos;
	CArrMulti<real, MAX_CELL_NUM, 2>ca2p;
	CArrMulti<real, MAX_CELL_NUM, 2>IP3;
public:

	CArr<CELL_STATE> state;

	CArr<CellIndex> fix_origin;
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

	CValue<real> zmax; //one value
	CValue<int> ncell;
	CValue<int> nder;
	CValue<int> nmemb;
	CValue<int> sw;
	CValue<int> next_uid;

	switcher swt;
	CellDataMan();
	~CellDataMan();
	__host__ thrust::device_ptr<CellPos> current_pos_arr();
	__host__ thrust::device_ptr<CellPos> next_pos_arr();
};

