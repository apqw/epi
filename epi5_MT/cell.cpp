#include "cell.h"
#include "cellmanager.h"
#include <iostream>
#include <istream>
#include <fstream>
#include <sstream>
#include <cinttypes>
#include <tbb/scalable_allocator.h>



double Cell::get_div_age_thresh(CELL_STATE state) {
	return state == FIX ? agki_max_fix
		: state == MUSUME ? agki_max
		: 0;
}

int Cell::correct_div_times(CELL_STATE state, int given_times) {
	return state == FIX ? div_max : given_times;
}



void Cell::set_index(size_t i)
{
	index = i;
}

void Cell::set_dermis(const Cell *const d)
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
	ca2p._migrate(_index);
	IP3._migrate(_index);
}


Cell::Cell(ctor_cookie,CELL_STATE _state, 
	SwapData<double[cont::MAX_CELL_NUM]>&ca2p_s,
	SwapData<double[cont::MAX_CELL_NUM]>&IP3_s,
	double _ex_inert,
	double _agek , double _ageb , double _ex_fat , double _in_fat, double _spr_nat_len,
	double _x, double _y, double _z,int _rest_div_time,
	double _radius , double _ca2p_avg ,
	bool _is_malignant ) :
	state(_state), 
	x(_x),
	y(_y),
	z(_z),
	ca2p(ca2p_s),
    ca2p_avg(_ca2p_avg),
	IP3(IP3_s),
	ex_inert(_ex_inert),
	agek(_agek),
	ageb(_ageb),
	ex_fat(_ex_fat),
	in_fat(_in_fat),
	spr_nat_len(_spr_nat_len),
    radius(_radius), 
	div_age_thresh(Cell::get_div_age_thresh(_state)),
	rest_div_times(Cell::correct_div_times(_state,_rest_div_time)),
	is_malignant(_is_malignant), diff_u(0) {
}

void* Cell::operator new(size_t s){
    return scalable_malloc(s);
}

void Cell::operator delete(void* p){
    scalable_free(p);
}
