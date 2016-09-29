#include "cell.h"

/** 細胞のステートに応じた分裂開始閾値を取得する */
real Cell::get_div_age_thresh(CELL_STATE state) {
    return state == FIX ? pm->agki_max_fix
        : state == MUSUME ? pm->agki_max
        : 0;
}

/** 細胞のステートに応じた分裂可能回数を計算 */
int Cell::correct_div_times(CELL_STATE state, int given_times) {
    return state == FIX ? pm->div_max : given_times;
}
/** dermisへの参照をセット */
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

/** dermisへの参照を取得 */
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


Cell::~Cell()
{
}
