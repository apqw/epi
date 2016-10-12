#include "cell.h"
#include <sstream>
/** �זE�̃X�e�[�g�ɉ���������J�n臒l���擾���� */
real Cell::get_div_age_thresh(CELL_STATE state) {
    return state == FIX ? pm->agki_max_fix
        : state == MUSUME ? pm->agki_max
        : 0;
}

/** �זE�̃X�e�[�g�ɉ���������\�񐔂��v�Z */
int Cell::correct_div_times(CELL_STATE state, int given_times) {
    return state == FIX ? pm->div_max : given_times;
}
/** dermis�ւ̎Q�Ƃ��Z�b�g */
void Cell::set_dermis(const Cell *const d)
{
    _dermis = d;
}

void Cell::set_index(size_t idx) {
    index = idx;
}

size_t Cell::get_index()const {
    return index;
}
void Cell::migrate(size_t destidx) {
    set_index(destidx);
}

/** dermis�ւ̎Q�Ƃ��擾 */
const Cell * Cell::dermis() const
{
    return _dermis;
}
Cell::Cell(ctor_cookie,CELL_STATE _state, int _fix_origin,
    real _ca2p,
    real _IP3,
    real _ex_inert,
    real _agek, real _ageb, real _ex_fat, real _in_fat, real _spr_nat_len ,
    real _x , real _y , real _z , int _rest_div_time ,
    real _radius, real _ca2p_avg,
    bool _is_malignant) :connected_cell(pm->MAX_CONNECT_CELL_NUM), state(_state),
    fix_origin((_state == DER || _state == MEMB) ? -1 : _fix_origin),
    x(_x),
    y(_y),
    z(_z),
    ca2p(_ca2p),
    ca2p_avg(_ca2p_avg),
    IP3(_IP3),
    ex_inert(_ex_inert),
    agek(_agek),
    ageb(_ageb),
    ex_fat(_ex_fat),
    in_fat(_in_fat),
    spr_nat_len(_spr_nat_len),
    radius(_radius),
    div_age_thresh(Cell::get_div_age_thresh(_state)),
    rest_div_times(Cell::correct_div_times(_state, _rest_div_time)),
    is_malignant(_is_malignant), diff_u(0)
{
}
std::string Cell::cell_info_str(){
	std::stringstream ss;
	ss<<"Cell information"<<std::endl
			<<"index:"<<index<<std::endl
			<<"state:"<<state_to_str(state)<<std::endl
			<<"x:"<<x()<<std::endl
			<<"y:"<<y()<<std::endl
			<<"z:"<<z()<<std::endl
			<<"radius:"<<radius<<std::endl
			<<"ca2p_avg:"<<ca2p_avg<<std::endl
			<<"agek:"<<agek<<std::endl
			<<"ageb:"<<ageb<<std::endl
			<<"ex_fat:"<<ex_fat<<std::endl
			<<"in_fat:"<<in_fat<<std::endl
			<<"spring nat len:"<<spr_nat_len<<std::endl
			<<"rest div times:"<<rest_div_times<<std::endl
			<<"fix_origin:"<<fix_origin<<std::endl;
	return ss.str();
}
Vec<3, real> Cell::pos_as_vector()const {
    return Vec<3, real>( x(),y(),z() );
}

