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

/*

	ManagerÇæÇ©ÇÁÇ∆åæÇ¡ÇƒÇ»ÇÒÇ≈Ç‡Ç©ÇÒÇ≈Ç‡ãlÇﬂçûÇ‹Ç»Ç¢Ç±Ç∆ÅIÅI

*/
class CellManager :Lockfree_push_stack<CellPtr, cont::MAX_CELL_NUM> {

    //std::vector<std::shared_ptr<Cell>> cell_store;

	SwapData<double[cont::MAX_CELL_NUM]>ca2p_s;
	SwapData<double[cont::MAX_CELL_NUM]>IP3_s;

	Lockfree_push_stack<size_t, cont::MAX_CELL_NUM> remove_queue;
	size_t nmemb, nder;
	size_t register_cell(const CellPtr& c);
	void _memb_init();
	void _load_from_file(std::string path);
    void _load_from_file_bin(std::string path);
	void init_internal(std::string init_data_path);
	std::atomic<uint_fast8_t> sw;
public:
	//void pos_swap();
	//void pos_copy(CellManager& cman);
	friend void ca2p_swap(CellManager& cman);
	friend void IP3_swap(CellManager& cman);
    friend void cman_init(CellManager& cells,const std::string& init_data_path,
                          bool use_last,const std::string& ld_uvp,const std::string& ld_w);
    friend void load_w(CellManager& cman,const std::string& ld_w);
    friend void load_uvp(CellManager& cman,const std::string& ld_uvp);
	friend void cell_pos_periodic_fix(CellManager& cman);
	friend void cornificate(CellManager& cman, Cell*const RESTRICT c);
	
    void output(const std::string& filename,bool binary_mode=false);
    void clean_up();

    CellPtr create(CELL_STATE _state, int stem_orig_id,double _x = 0, double _y = 0, double _z = 0,
		double _radius = cont::R_max, double _ca2p = cont::ca2p_init, double _ca2p_avg = cont::ca2p_init,
		double _IP3 = cont::IP3_init, double _ex_inert = cont::ex_inert_init,
		double _agek = 0, double _ageb = 0,
		double _ex_fat = 0, double _in_fat = 0,
		double _spr_nat_len = 0,
		int _rest_div_times = 0,
		bool _is_malignant = false);

	
	void add_remove_queue(size_t idx);
	void remove_exec();
	D_CELL_LOOP_ACCESSOR(all, 0, size());
	D_CELL_LOOP_ACCESSOR(memb, 0, nmemb);
	D_CELL_LOOP_ACCESSOR(der, nmemb, nder+nmemb);
	D_CELL_LOOP_ACCESSOR(other, nder + nmemb, size());
	D_CELL_LOOP_ACCESSOR(non_memb, nmemb, size());

	inline size_t size()const {
		return Lockfree_push_stack::size();
	}

	inline bool should_calc_ca() {
		return (sw >= cont::SW_THRESH);
	}

	inline void ca_calc_condition_reset() {
		sw = 0;
	}
};

void cman_init(CellManager& cells,const std::string& init_data_path,
               bool use_last=false,const std::string& ld_uvp="",const std::string& ld_w="");
void cell_pos_periodic_fix(CellManager& cman);
void ca2p_swap(CellManager& cman);
void IP3_swap(CellManager& cman);
void cornificate(CellManager& cman, Cell*const RESTRICT alive_cell);
void load_w(CellManager& cman,const std::string& ld_w);
void load_uvp(CellManager& cman,const std::string& ld_uvp);
void pos_copy(CellManager& cman);
