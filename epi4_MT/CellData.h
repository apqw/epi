#pragma once
#include "define.h"
#include "DVStore.h"
#include <unordered_map>

class CellData
{
	unsigned int _cell_num;
public:
	
	DVStore<double[MAX_CELL_NUM][3]> pos;
	DVStore<double[MAX_CELL_NUM]> radius;
	DVStore<double[MAX_CELL_NUM]> agek,ageb;
	DVStore<double[MAX_CELL_NUM]> ca2p;
	DVStore<double[MAX_CELL_NUM]> ca2p_avg;
	DVStore<double[MAX_CELL_NUM]> IP3;
	DVStore<double[MAX_CELL_NUM]> ex_inert;
	DVStore<double[MAX_CELL_NUM]> ex_fat,in_fat;
	DVStore<double[MAX_CELL_NUM]> nat_spring_len;
	double div_age_thresh[MAX_CELL_NUM];
	int rest_div_times[MAX_CELL_NUM];
	bool is_malignant[MAX_CELL_NUM];
	bool is_touch[MAX_CELL_NUM];
	int lat[MAX_CELL_NUM][3];
	double diffu[MAX_CELL_NUM];
	unsigned int connected_num[MAX_CELL_NUM];
	unsigned int connected_index[MAX_CELL_NUM][MAX_CONNECT_CELL_NUM];
	DVStore<std::unordered_map<unsigned int, double>[MAX_CELL_NUM]> gj;
	CELL_STATE state[MAX_CELL_NUM];

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

	void remove(unsigned int cell_idx);
	const unsigned int cell_num() const;

	CellData();
	~CellData();
};

