#include "CellManager.h"

#include <fstream>
#include <cinttypes>
#include <iostream>
std::vector<Cell*> CellLoadProc::tmp_pair;
CellManager::CellManager(size_t N):Lockfree_push_stack_dyn<Cell*>(N),remove_queue(N)
{
}

Cell* CellManager::create(CELL_STATE _state, int stem_orig_id, double _x, double _y, double _z, double _radius, double _ca2p, double _ca2p_avg, double _IP3, double _ex_inert, double _agek, double _ageb, double _ex_fat, double _in_fat, double _spr_nat_len, int _rest_div_times, bool _is_malignant)
{
    Cell* cptr=new Cell(Cell::ctor_cookie(),
        _state,
        stem_orig_id,
        _ca2p,
        _IP3,
        _ex_inert,
        _agek, _ageb, _ex_fat, _in_fat, _spr_nat_len,
        _x, _y, _z, _rest_div_times,

        _radius, _ca2p_avg, _is_malignant);
    cptr->set_index(this->push_back_with_index(cptr));
    return cptr;
}
void CellManager::add_remove_queue(size_t idx)
{
    remove_queue.push_back(idx);
}

void CellManager::remove_exec()
{
    for (size_t i = 0; i < remove_queue.size(); ++i) {
        size_t remove_idx = remove_queue[i];

        delete _data[remove_idx];
        _data[remove_idx] = _data[--_next];
        _data[remove_idx]->migrate(remove_idx);
    }
    remove_queue.clear();
}

CellLoadProc::CellLoadProc() {
    tmp_pair.resize(pm->MEMB_NUM_X*pm->MEMB_NUM_Y + 100000, nullptr);
}

void CellLoadProc::operator()(CellManager&cman, const CellTempStruct& cts) {
    auto cptr = cman.create(
        cts.state,
        cts.stem_orig_id,
        cts.x, cts.y, cts.z,
        cts.rad,
        cts.ca2p, cts.ca2p_avg,
        pm->IP3_init,
        pm->ex_inert_init,
        cts.agek, cts.ageb,
        cts.ex_fat, cts.fat,
        cts.spr_len,
        cts.div_times,
        (unsigned int)cts.stem_orig_id < pm->MALIG_NUM
    );

    if (tmp_pair.size() <= cts.id_count)tmp_pair.resize(cts.id_count + 1);
    if (cts.pair_cell_id > -1) {
        if (tmp_pair[cts.id_count] != nullptr) {

            if (tmp_pair[cts.id_count]->pair != nullptr) {
                throw std::logic_error("Wrong cell pairing");
            }
            cptr->pair = tmp_pair[cts.id_count];
            tmp_pair[cts.id_count]->pair = cptr;
        }
        else {
            tmp_pair.resize(cts.pair_cell_id + 1);
            tmp_pair[cts.pair_cell_id] = cptr;
        }
    }
}