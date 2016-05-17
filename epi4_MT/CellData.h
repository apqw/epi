#pragma once
#include "define.h"
#include "DVStore.h"
#include "atomic_double.h"
#include "LFStorage.h"
#include <tbb/parallel_for.h>
#include <tbb/parallel_for_each.h>
#include <tbb/blocked_range.h>
#include <unordered_map>
#include <atomic>


class CellData
{
	std::atomic_uint _cell_num;
	void remove_cell(unsigned int cell_idx);
	LFStorage<unsigned int, cont::MAX_CELL_NUM> remove_queue;
public:
	
	atomic_double pos[cont::MAX_CELL_NUM][3];
	double radius[cont::MAX_CELL_NUM];
	double agek[cont::MAX_CELL_NUM],ageb[cont::MAX_CELL_NUM];
	double ca2p[cont::MAX_CELL_NUM];
	double ca2p_avg[cont::MAX_CELL_NUM];
	double IP3[cont::MAX_CELL_NUM];
	double ex_inert[cont::MAX_CELL_NUM];
	double ex_fat[cont::MAX_CELL_NUM],in_fat[cont::MAX_CELL_NUM];
	double nat_spring_len[cont::MAX_CELL_NUM];
	double div_age_thresh[cont::MAX_CELL_NUM];
	int rest_div_times[cont::MAX_CELL_NUM];
	bool is_malignant[cont::MAX_CELL_NUM];
	bool is_touch[cont::MAX_CELL_NUM];
	int lat[cont::MAX_CELL_NUM][3];
	double diffu[cont::MAX_CELL_NUM];
	unsigned int connected_num[cont::MAX_CELL_NUM];
	unsigned int connected_index[cont::MAX_CELL_NUM][cont::MAX_CONNECT_CELL_NUM];
	std::unordered_map<unsigned int, double> gj[cont::MAX_CELL_NUM];
	CELL_STATE state[cont::MAX_CELL_NUM];

	//thread-safe
	unsigned int add(CELL_STATE _state,
		double _x = 0,
		double _y = 0,
		double _z = 0,
		double _radius = 0,
		double _agek = 0,
		double _ageb = 0,
		double _ca2p = 0,
		double _ca2p_avg = 0,
		double _IP3=0,
		double _ex_inert=0,
		double _ex_fat=0,
		double _in_fat=0,
		double _nat_spring_len=0,
		double _div_age_thresh=0,
		int _rest_div_times=0,
		bool _is_malignant=false,
		bool _is_touch=false);

	//thread-safe
	void queued_remove(unsigned int cell_idx);

	//non-thread-safe
	void exec_remove_queue();
	const unsigned int cell_num() const;
	
	template<class Fn>
	void foreach_parallel(Fn&& lmbd) {

		tbb::parallel_for(tbb::blocked_range<size_t>(0, _cell_num), [&](const tbb::blocked_range< size_t >& range) {
			for (size_t i = range.begin(); i != range.end(); ++i) {
				lmbd(i);
			}
		});
	}

	CellData();
	~CellData();
};

