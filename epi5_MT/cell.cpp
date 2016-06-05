#include "cell.h"
#include "cellmanager.h"
#include <iostream>
#include <istream>
#include <fstream>
#include <sstream>
#include <cinttypes>






size_t Cell::get_index() const
{
	return index;
}

void Cell::set_index(size_t i)
{
	index = i;
}

void Cell::set_dermis(Cell * d)
{
	_dermis = d;
}

const Cell * Cell::dermis() const
{
	return _dermis;
}

void Cell::migrate(size_t dest_idx)
{
	set_index(dest_idx);
	size_t _index = get_index();
	x._migrate(_index);
	y._migrate(_index);
	z._migrate(_index);
	ca2p._migrate(_index);
	IP3._migrate(_index);
	ex_inert._migrate(_index);
	agek._migrate(_index);
	ageb._migrate(_index);
	in_fat._migrate(_index);
	ex_fat._migrate(_index);
	spr_nat_len._migrate(_index);
	gj._migrate(_index);
	//rest_div_times._migrate(_index);
}


Cell::Cell(ctor_cookie,CELL_STATE _state, 
	SwapData<atomic_double[cont::MAX_CELL_NUM][3]>& pos_s,
	SwapData<double[cont::MAX_CELL_NUM]>&ca2p_s,
	SwapData<double[cont::MAX_CELL_NUM]>&IP3_s,
	SwapData<double[cont::MAX_CELL_NUM]>&ex_inert_s,
	SwapData<double[cont::MAX_CELL_NUM]>&agek_s,
	SwapData<double[cont::MAX_CELL_NUM]>&ageb_s,
	SwapData<double[cont::MAX_CELL_NUM]>&ex_fat_s,
	SwapData<double[cont::MAX_CELL_NUM]>&in_fat_s,
	SwapData<double[cont::MAX_CELL_NUM]>&spr_nat_len_s,
	SwapData<std::unordered_map<Cell*, double>[cont::MAX_CELL_NUM]>&gj_s,
	double _radius , double _ca2p_avg ,
	double _div_age_thresh ,
	bool _is_malignant ) :
	state(_state), 
	x(pos_s),
	y(pos_s),
	z(pos_s),
	ca2p(ca2p_s),
	IP3(IP3_s),
	ex_inert(ex_inert_s),
	agek(agek_s),
	ageb(ageb_s),
	ex_fat(ex_fat_s),
	in_fat(in_fat_s),
	spr_nat_len(spr_nat_len_s),
	gj(gj_s),
	radius(_radius), ca2p_avg(_ca2p_avg), div_age_thresh(_div_age_thresh), is_malignant(_is_malignant), diff_u(0) {
}
