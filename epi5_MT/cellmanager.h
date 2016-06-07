#pragma once
#include <memory>
#include <unordered_map>
#include "atomics.h"
#include "define.h"
#include "swapdata.h"
#include <tbb/parallel_for.h>
#include <tbb/parallel_for_each.h>
#include <tbb/blocked_range.h>
#include <vector>


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
	tbb::parallel_for_each(_data+(start),_data+ (end), lmbd);\
}


class CellManager :public Lockfree_push_stack<CellPtr, cont::MAX_CELL_NUM> {

	std::vector<std::shared_ptr<Cell>> cell_store;

	SwapData<double[cont::MAX_CELL_NUM]>ca2p_s;
	SwapData<double[cont::MAX_CELL_NUM]>IP3_s;
	SwapData<double[cont::MAX_CELL_NUM]>ex_inert_s;

	Lockfree_push_stack<size_t, cont::MAX_CELL_NUM> remove_queue;
	size_t nmemb, nder;
public:
	void pos_swap();
	void pos_copy();
	void ca2p_swap();
	void IP3_swap();
	void ex_inert_swap();
	void agek_swap();
	void ageb_swap();
	void ex_fat_swap();
	void in_fat_swap();
	void spr_nat_len_swap();
	friend void cman_load_from_file(CellManager& cman, std::string path);
	friend void cell_pos_periodic_fix(CellManager& cman);
	

	CellPtr create(CELL_STATE _state, double _x = 0, double _y = 0, double _z = 0,
		double _radius = cont::R_max, double _ca2p = cont::ca2p_init, double _ca2p_avg = cont::ca2p_init,
		double _IP3 = cont::IP3_init, double _ex_inert = cont::ex_inert_init,
		double _agek = 0, double _ageb = 0,
		double _ex_fat = 0, double _in_fat = 0,
		double _spr_nat_len = 0,
		double _div_age_thresh = 0,
		int _rest_div_times = 0,
		bool _is_malignant = false);

	size_t register_cell(const CellPtr& c);
	void add_remove_queue(size_t idx);
	void remove_exec();
	D_CELL_LOOP_ACCESSOR(all, 0, size());
	D_CELL_LOOP_ACCESSOR(memb, 0, nmemb);
	D_CELL_LOOP_ACCESSOR(der, nmemb, nder+nmemb);
	D_CELL_LOOP_ACCESSOR(other, nder + nmemb, size());
	D_CELL_LOOP_ACCESSOR(non_memb, nmemb, size());

	std::atomic<uint_fast8_t> sw;
};

void cman_load_from_file(CellManager& cman, std::string path);
void cell_pos_periodic_fix(CellManager& cman);
