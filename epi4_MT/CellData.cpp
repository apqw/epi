#include "CellData.h"

#include <cassert>
//thread-safe
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
	assert(_cell_num < cont::MAX_CELL_NUM);
	unsigned int my_cell_num = _cell_num++;
	state[my_cell_num] = _state;
	pos[my_cell_num][0] = _x;
	pos[my_cell_num][1] = _y;
	pos[my_cell_num][2] = _z;
	radius[my_cell_num] = _radius;
	agek[my_cell_num] = _agek;
	ageb[my_cell_num] = _ageb;
	ca2p[my_cell_num] = _ca2p;
	ca2p_avg[my_cell_num] = _ca2p_avg;
	IP3[my_cell_num] = _IP3;
	ex_inert[my_cell_num] = _ex_inert;
	ex_fat[my_cell_num] = _ex_fat;
	in_fat[my_cell_num] = _in_fat;
	nat_spring_len[my_cell_num] = _nat_spring_len;
	div_age_thresh[my_cell_num] = _div_age_thresh;
	rest_div_times[my_cell_num] = _rest_div_times;
	is_malignant[my_cell_num] = _is_malignant;
	is_touch[my_cell_num] = _is_touch;
	gj[my_cell_num].clear();
	return my_cell_num;
}
//not thread safe
void CellData::remove_cell(unsigned int cidx)
{
	assert(_cell_num > 0 && cidx<_cell_num);
	_cell_num--;
	state[cidx] = state[_cell_num];
	pos[cidx][0] = (double)pos[_cell_num][0];
	pos[cidx][1] = (double)pos[_cell_num][1];
	pos[cidx][2] = (double)pos[_cell_num][2];
	radius[cidx] = radius[_cell_num];
	agek[cidx] = agek[_cell_num];
	ageb[cidx] = ageb[_cell_num];
	ca2p[cidx] = ca2p[_cell_num];
	ca2p_avg[cidx] = ca2p_avg[_cell_num];
	IP3[cidx] = IP3[_cell_num];
	ex_inert[cidx] = ex_inert[_cell_num];
	ex_fat[cidx] = ex_fat[_cell_num];
	in_fat[cidx] = in_fat[_cell_num];
	nat_spring_len[cidx] = nat_spring_len[_cell_num];
	div_age_thresh[cidx] = div_age_thresh[_cell_num];
	rest_div_times[cidx] = rest_div_times[_cell_num];
	is_malignant[cidx] = is_malignant[_cell_num];
	is_touch[cidx] = is_touch[_cell_num];
	gj[cidx] = gj[_cell_num];
	for (int i = 0; i < _cell_num; i++) {
		if (gj[i].count(_cell_num)>0) {
			gj[i][cidx] = gj[i].at(_cell_num);
		}
	}
}

void CellData::queued_remove(unsigned int cell_idx)
{
	remove_queue.push_back(cell_idx);
}

void CellData::exec_remove_queue()
{
	unsigned int sz = remove_queue.count();
	for (unsigned int i = 0; i < sz; i++) {
		remove_cell(remove_queue.data()[i]);
	}
	remove_queue.clear();
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
