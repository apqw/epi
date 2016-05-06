#pragma once
#include "define.h"
#include "component.h"
#include "Cell.h"
#include <istream>
#include <fstream>
class Field
{
private:

	const int _MAX_CELL_NUM = 30000;

	CellMan cells;
	Arr3D<DV<double>> ATP;
	Arr3D<DV<double>> ext_stim;
	Arr3D<Cell*> cell_map;
	Arr3D<int> cell_map2;
	Arr3D<int> air_stim_flg;
	Arr3D<double> age_map;
	Arr3D<int> cornif_map;
	void interact_cell();
	void cell_state_renew();
	void cell_pos_periodic_fix();
	void connect_cells();
public:
	Field(
		int __MAX_CELL_NUM=30000
		):_MAX_CELL_NUM(__MAX_CELL_NUM),
		ATP(Arr3D<DV<double>>(cont::NX+1, Arr2D<DV<double>>(cont::NY+1, Arr1D<DV<double>>(cont::NZ+1, cont::ATP_init)))),
		ext_stim(Arr3D<DV<double>>(cont::NX+1, Arr2D<DV<double>>(cont::NY+1, Arr1D<DV<double>>(cont::NZ+1, cont::ext_stim_init)))),
		cell_map(Arr3D<Cell*>(cont::NX + 1, Arr2D<Cell*>(cont::NY + 1, Arr1D<Cell*>(cont::NZ + 1, nullptr)))),
		cell_map2(Arr3D<int>(cont::NX + 1, Arr2D<int>(cont::NY + 1, Arr1D<int>(cont::NZ + 1, 0)))),
		air_stim_flg(Arr3D<int>(cont::NX + 1, Arr2D<int>(cont::NY + 1, Arr1D<int>(cont::NZ + 1, 0)))),
		age_map(Arr3D<double>(cont::NX + 1, Arr2D<double>(cont::NY + 1, Arr1D<double>(cont::NZ + 1, 0)))),
		cornif_map(Arr3D<int>(cont::NX + 1, Arr2D<int>(cont::NY + 1, Arr1D<int>(cont::NZ + 1, 0))))
		{
	}


	void init_with_file(std::ifstream& dstrm);
	void cell_dynamics();
	void main_loop();
	void setup_map();
	void calc_b();
	
	//~Field();
};

