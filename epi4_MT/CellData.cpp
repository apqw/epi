#include "CellData.h"

#include <cassert>

unsigned int CellData::add(CELL_STATE _state,
	double _x = 0,
	double _y = 0,
	double _z = 0,
	double _radius = 0,
	double _agek = 0,
	double _ageb = 0,
	double _ca2p = 0,
	double _ca2p_avg = 0,
	double _IP3 = 0,
	double _ex_inert = 0,
	double _ex_fat = 0,
	double _in_fat = 0,
	double _nat_spring_len = 0,
	double _div_age_thresh = 0,
	int _rest_div_times = 0,
	bool _is_malignant = false,
	bool _is_touch = false){
	assert(_cell_num < MAX_CELL_NUM);
	state[_cell_num] = _state;
	pos._raw_old()[_cell_num][0] = _x;
	pos._raw_old()[_cell_num][1] = _y;
	pos._raw_old()[_cell_num][2] = _z;
	radius._raw_old()[_cell_num] = _radius;
	agek._raw_old()[_cell_num] = _agek;
	ageb._raw_old()[_cell_num] = _ageb;
	ca2p._raw_old()[_cell_num] = _ca2p;
	ca2p_avg._raw_old()[_cell_num] = _ca2p_avg;
	IP3._raw_old()[_cell_num] = _IP3;
	ex_inert._raw_old()[_cell_num] = _ex_inert;
	ex_fat._raw_old()[_cell_num] = _ex_fat;
	in_fat._raw_old()[_cell_num] = _in_fat;
	nat_spring_len._raw_old()[_cell_num] = _nat_spring_len;
	div_age_thresh[_cell_num] = _div_age_thresh;
	rest_div_times[_cell_num] = _rest_div_times;
	is_malignant[_cell_num] = _is_malignant;
	is_touch[_cell_num] = _is_touch;
	gj._raw_old()[_cell_num].clear();
	return _cell_num++;
}

void CellData::remove(unsigned int cidx)
{
	assert(_cell_num > 0 && cidx<_cell_num);
	--_cell_num;
	state[cidx] = state[_cell_num];
	pos._raw_old()[cidx][0] = pos._raw_old()[_cell_num][0];
	pos._raw_old()[cidx][1] = pos._raw_old()[_cell_num][1];
	pos._raw_old()[cidx][2] = pos._raw_old()[_cell_num][2];
	radius._raw_old()[cidx] = radius._raw_old()[_cell_num];
	agek._raw_old()[cidx] = agek._raw_old()[_cell_num];
	ageb._raw_old()[cidx] = ageb._raw_old()[_cell_num];
	ca2p._raw_old()[cidx] = ca2p._raw_old()[_cell_num];
	ca2p_avg._raw_old()[cidx] = ca2p_avg._raw_old()[_cell_num];
	IP3._raw_old()[cidx] = IP3._raw_old()[_cell_num];
	ex_inert._raw_old()[cidx] = ex_inert._raw_old()[_cell_num];
	ex_fat._raw_old()[cidx] = ex_fat._raw_old()[_cell_num];
	in_fat._raw_old()[cidx] = in_fat._raw_old()[_cell_num];
	nat_spring_len._raw_old()[cidx] = nat_spring_len._raw_old()[_cell_num];
	div_age_thresh[cidx] = div_age_thresh[_cell_num];
	rest_div_times[cidx] = rest_div_times[_cell_num];
	is_malignant[cidx] = is_malignant[_cell_num];
	is_touch[cidx] = is_touch[_cell_num];
	gj._raw_old()[cidx] = gj._raw_old()[_cell_num];
	for (int i = 0; i < _cell_num; i++) {
		if (gj._raw_old()[i].count(_cell_num)>0) {
			gj._raw_old()[i][cidx] = gj._raw_old()[i].at(_cell_num);
		}
	}
}

const unsigned int CellData::cell_num() const
{
	return _cell_num;
}

CellData::CellData() :pos()
{
}


CellData::~CellData()
{
}
