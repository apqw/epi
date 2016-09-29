#pragma once
#include "cell.h"
#include <tbb/parallel_for.h>
#include <tbb/parallel_for_each.h>
#include <tbb/blocked_range.h>

#define D_CELL_LOOP_ACCESSOR(prefix,start,end)\
template<class Fn>\
inline void prefix##_foreach_with_index(const Fn& lmbd) {\
	for (size_t i = (start); i < (end); ++i) {\
		lmbd((*this)[i], i);\
	}\
}\
\
template<class Fn>\
inline void prefix##_foreach(const Fn& lmbd) {\
	for (size_t i = (start); i < (end); ++i) {\
		lmbd((*this)[i]);\
	}\
}\
\
template<class Fn>\
inline void prefix##_foreach_parallel(const Fn& lmbd) {\
	tbb::parallel_for<size_t,Fn>((start), (end), lmbd);\
}\
template<class Fn>\
inline void prefix##_foreach_parallel_native(const Fn& lmbd) {\
	\
	tbb::parallel_for_each(&_data[0]+(start),&_data[0]+ (end), lmbd);\
}

class CellManager:public Lockfree_push_stack_dyn<Cell*>
{
    Lockfree_push_stack_dyn<size_t> remove_queue;
    size_t nmemb, nder;
    void memb_init();
public:
    CellManager(size_t);
    Cell* create(CELL_STATE _state, int stem_orig_id, real _x = 0, real _y = 0, real _z = 0,
        real _radius = pm->R_max, real _ca2p = pm->ca2p_init, real _ca2p_avg = pm->ca2p_init,
        real _IP3 = pm->IP3_init, real _ex_inert = pm->ex_inert_init,
        real _agek = 0, real _ageb = 0,
        real _ex_fat = 0, real _in_fat = 0,
        real _spr_nat_len = 0,
        int _rest_div_times = 0,
        bool _is_malignant = false);
    void add_remove_queue(size_t idx);
    void remove_exec();
    void load(const std::string&);
    D_CELL_LOOP_ACCESSOR(all, 0, size());
    D_CELL_LOOP_ACCESSOR(memb, 0, nmemb);
    D_CELL_LOOP_ACCESSOR(der, nmemb, nder + nmemb);
    D_CELL_LOOP_ACCESSOR(other, nder + nmemb, size());
    D_CELL_LOOP_ACCESSOR(non_memb, nmemb, size());

    //~CellManager();
};

