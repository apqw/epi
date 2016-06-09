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
#include "ca2p.h"
#include <tbb/task_group.h>

inline void cell_dynamics(CellManager & cellset) {
	cell_interaction(cellset);
	pos_copy(cellset);
	cell_state_renew(cellset);
	
	//does copy
	cell_pos_periodic_fix(cellset);
	

	connect_cell(cellset);
}

double calc_zzmax(CellManager& cman) {
	double zmax = 0;
	cman.other_foreach([&zmax](Cell*const RESTRICT c) {
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
	init_precalc_lat();
	init_precalc_per();

	bool flg_forced_sc = FORCE_CORNIF;
	int num_sc = 0;
	auto ATP=std::make_unique<SwapData<FArr3D<double>>>();
	auto ext_stim= std::make_unique<SwapData<FArr3D<double>>>();
	auto cell_map1 = std::make_unique<FArr3D<const Cell*>>();
	auto cell_map2 = std::make_unique<FArr3D<uint_fast8_t>>();
	double zzmax = 0;
    printf("current cell num:%zd\n", cellset->size());
	auto& cman = *cellset;
    for (size_t i = 0; i < NUM_ITR; i++) {
        if(i%100==0)printf("loop:%zd\n", i);
			cell_dynamics(cman);
			zzmax = calc_zzmax(cman);
			setup_map_lat(cman, *cell_map1, *cell_map2);
		
		calc_ext_stim(*ext_stim, *cell_map1, *cell_map2, zzmax);
		if (i*DT_Cell > T_TURNOVER&&flg_forced_sc) {
			flg_forced_sc = false;
			printf("forced cornif\n");
			initialize_sc(cman,zzmax);
			num_sc = NUM_SC_INIT;
		}

		if (cman.should_calc_ca() || num_sc > 0) {
			printf("calc ca...\n");
			calc_ca2p(cman, *ATP, ext_stim->first(), *cell_map1, *cell_map2, zzmax);
			//cells.all_cell_update();
			printf("end.\n");
			if (num_sc>0)num_sc--;
			cman.ca_calc_condition_reset();
		}
	}
}
