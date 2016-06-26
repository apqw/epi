#include "cell.h"
#include "cellmanager.h"
#include <iostream>
#include <istream>
#include <fstream>
#include <sstream>
#include <cinttypes>
#include <tbb/scalable_allocator.h>

/** 細胞のステートに応じた分裂開始閾値を取得する */
double Cell::get_div_age_thresh(CELL_STATE state) {
	return state == FIX ? agki_max_fix
		: state == MUSUME ? agki_max
		: 0;
}

/** 細胞のステートに応じた分裂可能回数を計算 */
int Cell::correct_div_times(CELL_STATE state, int given_times) {
	return state == FIX ? div_max : given_times;
}


/** インデックスをセット */
void Cell::set_index(size_t i)
{
	index = i;
}

/** dermisへの参照をセット */
void Cell::set_dermis(const Cell *const d)
{
	_dermis = d;
}

/** dermisへの参照を取得 */
const Cell * Cell::dermis() const
{
	return _dermis;
}

/** 内部データを指定されたインデックスに応じて移動 */
void Cell::migrate(size_t dest_idx)
{
	set_index(dest_idx);
	size_t _index = get_index();
	ca2p._migrate(_index);
	IP3._migrate(_index);
}

/** ctor */
Cell::Cell(ctor_cookie,CELL_STATE _state, int _fix_origin,
	SwapData<double[cont::MAX_CELL_NUM]>&ca2p_s,
	SwapData<double[cont::MAX_CELL_NUM]>&IP3_s,
	double _ex_inert,
	double _agek , double _ageb , double _ex_fat , double _in_fat, double _spr_nat_len,
	double _x, double _y, double _z,int _rest_div_time,
	double _radius , double _ca2p_avg ,
	bool _is_malignant ) :
	state(_state), 
    fix_origin((_state==DER||_state==MEMB)?-1:_fix_origin),
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
