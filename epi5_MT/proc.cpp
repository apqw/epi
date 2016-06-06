#include "proc.h"
#include "swapdata.h"
#include "Field.h"
#include "cell.h"
#include "cellmanager.h"
#include "define.h"
#include "cell_init.h"
#include "cell_interaction.h"
#include "cell_connection.h"
#include "cell_state_renew.h"
#include "map.h"
#include "utils.h"
#include "ext_stim.h"
#include <tbb/task_group.h>

inline void cell_dynamics(CellManager & cellset) {
	cell_interaction(cellset);
	cell_state_renew(cellset);


	cell_pos_periodic_fix(cellset);
	cellset.ageb_swap();
	cellset.agek_swap();
	cellset.in_fat_swap();
	cellset.ex_fat_swap();
	cellset.spr_nat_len_swap();
	cellset.pos_copy();

	connect_cell(cellset);
}

double calc_zzmax(CellManager& cman) {
	double zmax = 0;
	cman.other_foreach([&](Cell*c) {
		auto& st = c->state;
		if (st==DEAD||st==ALIVE||st==MUSUME||st==FIX){//get_state_mask(st)&(DEAD_M | ALIVE_M | MUSUME_M | FIX_M)) {
			if (zmax < c->z())zmax = c->z();
		}
	});
	return zmax;
}

void proc(std::string init_data_path,std::string param_path,std::string init_uvp_data,std::string init_ATP_data,std::string init_ext_stim_data)
{
	using namespace cont;
	auto cellset = std::make_unique<CellManager>();
	cman_init(*cellset, init_data_path);
	init_cell_interaction(); //matomeru?
	init_cell_state_renewal();
	init_precalc_lat();
	init_precalc_per();

	auto ATP=std::make_unique<SwapData<Field<double,NX+1, NY+1, NZ+1>>>();
	auto ext_stim= std::make_unique<SwapData<Field<double,NX+1, NY+1, NZ+1>>>();
	auto cell_map1 = std::make_unique<Field<Cell*, NX+1, NY+1, NZ+1>>();
	auto cell_map2 = std::make_unique<Field<uint_fast8_t, NX + 1, NY + 1, NZ + 1>>();
	double zzmax = 0;
	printf("current cell num:%d\n", cellset->size());
	for (int i = 0; i < 10000; i++) {
		if(i%100==0)printf("loop:%d\n", i);
			cell_dynamics(*cellset);
			zzmax = calc_zzmax(*cellset);
			setup_map_lat(*cellset, *cell_map1, *cell_map2);
		
		calc_ext_stim(*ext_stim, *cell_map1, *cell_map2, zzmax);
	}
}
