#pragma once
#include <memory>
#include <array>
#include "define.h"
#include "component.h"
#include "util_func.h"
#include <unordered_map>
#include <tbb/parallel_for.h>
#include <tbb/parallel_for_each.h>
#include <tbb/blocked_range.h>

//class Cell;
using CellPtr = std::shared_ptr<Cell>;
/*
	IMPLEMENTS:
		ORDER:MEMB DER others...


*/

template<unsigned N>
class LWCellMan {
	std::array<Cell*, N> cell;
	std::atomic<int> _count;
public:

	LWCellMan() :_count(0) {};
	void add(const CellPtr& c) {
		cell[_count++] = c.get(); //_count++;
	}

	void add(Cell* cptr) {
		cell[_count++] = cptr; //_count++;
	}
	void reset() {
		count = 0;
	}

	void set_count(int c) {
		_count = c;
	}
	int count()const {
		return _count;
	}
	template<class L>
    void foreach(const L& lmbd) {
		for (int i = 0; i < _count; i++) {
			lmbd(cell[i]);
		}
	}

	bool exist(Cell *const  c) {
		for (int i = 0; i < _count; ++i) {
			if (cell[i] == c)return true;
		}
		return false;
	}

    std::array<Cell*, N>& _cell(){
        return cell;
    }
	

};

struct MEMB_bend_data {
	CellPtr memb_u; CellPtr memb_b; CellPtr memb_l; CellPtr memb_r;
	CellPtr memb_bb; CellPtr memb_ll;
	Vec3<double> nv; Vec3<double> mv;
	double ipn, ipm; double dn, dm;

	Vec3<double> memb_bend_force_sqr();
};

class Cell
{
private:

	std::array<bool,cont::STATE_NUM> f_interacted;
	template<class cf>
	void _interact(Cell* me, Cell* oppo) {
		double ljm = cf()(me, oppo);
		if (fabs(ljm) >= 100) {
			printf("error too strong interaction: ljm=%lf\n", ljm);
			assert(fabs(ljm) < 100);
		}
		
		auto dum = (cont::DT_Cell*ljm)*p_diff_v3(me->pos, oppo->pos);
	
		me->pos += dum;
		oppo->pos -= dum;
	}
	struct c_der_to_der {
		double operator()(Cell* me, Cell* oppo);
	};
	struct c_al_air_de_to_al_air_de_fix_mu {
		double operator()(Cell* me, Cell* oppo);
	};
	struct c_fix_mu_to_fix_mu {
		double operator()(Cell* me, Cell* oppo);
	};
	struct c_mu_to_memb {
		double operator()(Cell* me, Cell* oppo);
	};
	struct c_fix_to_memb {
		double operator()(Cell* me, Cell* oppo);
	};
	struct c_memb_to_memb {
		double operator()(Cell* me, Cell* oppo);
	};

	struct c_other {
		double operator()(Cell* me, Cell* oppo);
	};
	/*
	struct wall_interact {
		void operator()(Cell* me);
	};
	
	*/
	double ageb_const();
	double agek_const();
	double k_lipid_release();
	double k_lipid();
	bool pending_kill = false;
	bool pair_generated = false;
public:


	void DER_interact();

	void AL_AIR_DE_interact();

	void FIX_interact();

	void MUSUME_interact();

	void MEMB_interact();

	void wall_interact();

	void pair_interact();
	/*
		dirty impl
	*/
	void memb_bend_calc1();
	void memb_bend_calc2();

	void memb_bend_interact();

	void set_dermis();

	bool divide_try();

	void FIX_state_renew();
	void MUSUME_state_renew();
	void DEAD_AIR_state_renew();
	void ALIVE_state_renew();
	
	void set_as_deleted();
	void set_as_pair_gen();

	bool should_deleted();
	bool has_new_pair();
	void set_as_no_more_new_pair();

	void pair_disperse();
	void set_lattice();

	LWCellMan<400> connected_cell;
	lazy_update<CELL_STATE> state;
	CellPtr pair;
	Cell* dermis;
	Vec3<DV<double>> pos = { 0.0,0.0,0.0 };
	DV<double> radius = 0;
	DV<double> ca2p = 0;
	DV<double> ca2p_avg = 0;
	DV<double> IP3 = 0;
	DV<double> ex_inert = 0;
	DV<double> agek = 0, ageb = 0;
	DV<double> ex_fat = 0, in_fat = 0;
	DV<double> spring_nat_len = 0;
	DV<double> spring_force = 0;
	DV<double> div_age_thresh = 0;
	DV<double> poisson_div_thresh = 0;
	DV<int> rest_div_times = 0;
	DV<bool> is_malignant = 0;
    DV<bool> is_touch = 0;
	Vec3<int> lat;
	std::unordered_map<Cell*,DV<double>> gj;
	double diffu = 0; //tmp
    int __my_index=0;
	
	
	MEMB_bend_data mbd;

	Cell(
		CELL_STATE _state = UNUSED,
		Vec3<DV<double>> _pos = { 0,0,0 },
		DV<double> _radius = 0,
		DV<double> _ca2p = 0,
		DV<double> _ca2p_avg = 0,
		DV<double> _IP3 = 0,
		DV<double> _ex_inert = 0,
		DV<double> _agek = 0, DV<double>_ageb = 0,
		DV<double> _ex_fat = 0, DV<double> _in_fat = 0,
		DV<double> _spring_nat_len = 0, DV<double> _spring_force = 0,
		DV<double> _div_age_thresh = 0, DV<double> _poisson_div_thresh = 0,
		DV<int> _rest_div_times = 0, DV<bool> _is_malignant = false,
		DV<bool> _is_touch = false) :
		state(_state), pos(_pos), radius(_radius), ca2p(_ca2p), ca2p_avg(_ca2p_avg),
		IP3(_IP3), ex_inert(_ex_inert), agek(_agek), ageb(_ageb), ex_fat(_ex_fat),
		in_fat(_in_fat), spring_nat_len(_spring_nat_len),
		spring_force(_spring_force), div_age_thresh(_div_age_thresh),
		poisson_div_thresh(_poisson_div_thresh), rest_div_times(_rest_div_times),
		is_malignant(_is_malignant), is_touch(_is_touch) {

	}

	void do_interact_w_connected();
	void update();
	//Cell();
	~Cell();
};



class CellMan {
	std::vector<CellPtr> cell;
	std::vector<int> remove_index;
	std::vector<CellPtr> add_cell;
	int memb_num; int der_num;
	bool nmemb_is_set; int nder_is_set;
public:
	void set_memb_num(int n) {
		assert(!nmemb_is_set);
		memb_num = n;
		nmemb_is_set = true;
	}

	void set_der_num(int n) {
		assert(!nder_is_set);
		der_num = n;
		nder_is_set = true;
	}


    CellMan():nmemb_is_set(false),nder_is_set(false) {
		cell.reserve(cont::DEFAULT_MAX_CELL_NUM);
	};

	void add_direct(CellPtr& made_ptr) {
		cell.push_back(made_ptr);
	}
    void add_queue(const CellPtr& made_ptr) {
		add_cell.push_back(made_ptr);
	}

	void remove_queue(int idx) {
		remove_index.push_back(idx);
	}

	std::vector<CellPtr>& _raw_cell_set() {
		return cell;
	}



	/*
	remove instantly
	*/
	template<class L>
    void remove_if(const L& pred) {
		auto it = std::remove_if(cell.begin(), cell.end(), pred);
		cell.erase(it, cell.end());
	}

	template<class L>
    void remove_if_queue(const L& pred) {
		size_t sz = cell.size();
		for (int i = 0; i < sz; ++i) {
			if (pred(cell[i], i)) {
				remove_queue(i);
			}
		}
	}

	template<class L>
    void foreach(const L& lmbd) {
		size_t sz = cell.size();
		
		for (int i = 0; i < sz; ++i) {
			lmbd(cell[i], i);
		}
		
	}

	template<class L>
    void foreach_parallel(const L& lmbd) {
		size_t sz = cell.size();
		tbb::parallel_for(tbb::blocked_range<size_t>(0, sz), [&](const tbb::blocked_range< size_t >& range) {
			for (size_t i = range.begin(); i != range.end(); ++i) {
				lmbd(cell[i], i);
			}
		});
	}

	template<class L>
	void foreach_parallel_native(const L& lmbd) {
		tbb::parallel_for_each(cell.begin(), cell.end(), lmbd);
	}


	template<class L>
    void memb_foreach(const L& lmbd) {
		assert(nmemb_is_set);

		
		for (int i = 0; i < memb_num; ++i) {
		lmbd(cell[i], i);
		}
		
	}

	template<class L>
    void memb_foreach_parallel(const L& lmbd) {
		assert(nmemb_is_set);

		tbb::parallel_for(tbb::blocked_range<int>(0, memb_num), [&](const tbb::blocked_range< int >& range) {
			for (int i = range.begin(); i != range.end(); ++i) {
				lmbd(cell[i], i);
			}
		});
	}

	template<class L>
	void memb_foreach_parallel_native(const L& lmbd) {
		assert(nmemb_is_set);
		tbb::parallel_for_each(cell.begin(), cell.begin() + memb_num, lmbd);
		/*
		tbb::parallel_for(tbb::blocked_range<int>(0, memb_num), [&](const tbb::blocked_range< int >& range) {
			for (int i = range.begin(); i != range.end(); ++i) {
				lmbd(cell[i], i);
			}
		});
		*/
	}

	template<class L>
    void der_foreach(const L& lmbd) {
		assert(nder_is_set);
		for (int i = memb_num + 1; i < memb_num + der_num; ++i) {
			lmbd(cell[i], i);
		}
	}

	template<class L>
    void other_foreach(const L& lmbd) {
		assert(nmemb_is_set);
		assert(nder_is_set);
		size_t sz = cell.size();
		for (int i = memb_num +der_num; i < sz; ++i) {
			lmbd(cell[i], i);
		}
	}

	template<class L>
	void other_foreach_parallel(const L& lmbd) {
		assert(nmemb_is_set);
		assert(nder_is_set);
		size_t sz = cell.size();
		tbb::parallel_for(tbb::blocked_range<int>(memb_num + der_num , sz), [&](const tbb::blocked_range< int >& range) {
			for (int i = range.begin(); i != range.end(); ++i) {
				lmbd(cell[i], i);
			}
		});
	}

	template<class L>
	void other_foreach_parallel_native(const L& lmbd) {
		assert(nmemb_is_set);
		assert(nder_is_set);
		tbb::parallel_for_each(cell.begin() + memb_num + der_num, cell.end(), lmbd);
	}

	void all_cell_update() {
        foreach_parallel([](CellPtr& c, size_t i) {
			c->update();
		});
	}
	void update() {
		for (auto ridx : remove_index) {
			cell[ridx] = cell.back(); //delete(overwrite)
			cell.pop_back();
		}
		cell.insert(cell.end(), add_cell.begin(), add_cell.end());
		add_cell.clear();
		remove_index.clear();
	}
};
