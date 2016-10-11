#include "calc.h"
#include "../global.h"
#include "../define.h"
#include "cpu2/cell_connection.h"
#include "cpu2/CellManager.h"
#include "cpu2/cell_interaction.h"
#include "cpu2/cell_state_renew.h"
#include "../fsys.h"
#include "../utils.h"
#include "cpu2/map.h"
#include <cstdint>
static void output_cell_data(CellManager&cman,int index){
    std::cout << "output:" << index << std::endl;
	cman.output(pm->outputdir+"/"+std::to_string(index));
}

static void initialize(CellManager&cman){
    cman._memb_init();
	connect_cell(cman);
	cman.init_value();

	init_precalc_lat();
	init_precalc_per();

    make_dir(pm->outputdir.c_str());
}



static void cell_dynamics(CellManager&cman) {
    cell_interaction(cman);
    cman.pos_update();
    cell_state_renew(cman);
    cman.pos_periodic_fix();
    connect_cell(cman);
}

static real calc_zzmax(CellManager& cman) {
	real zmax = 0;
	cman.other_foreach([&zmax](Cell*const RESTRICT c) {
		auto& st = c->state;
		if (st==DEAD||st==ALIVE||st==MUSUME||st==FIX){//get_state_mask(st)&(DEAD_M | ALIVE_M | MUSUME_M | FIX_M)) {
			if (zmax < c->z())zmax = c->z();
		}
	});
	return zmax;
}

void calc_with_cpu(CellManager&cman){
    initialize(cman);

    real zzmax=0.0;
    auto cell_map1 = make_unique_c11<Dyn3DArr<const Cell*>>(pm->NX+1,pm->NY+1,pm->NZ+1);
    auto cell_map2 = make_unique_c11<Dyn3DArr<uint_fast8_t>>(pm->NX+1,pm->NY+1,pm->NZ+1);

	for(size_t i=0;i<pm->NUM_ITR;i++){
		if(i%pm->CUT==0){
			output_cell_data(cman, static_cast<int>(i/pm->CUT));
		}
        cell_dynamics(cman);
        zzmax=calc_zzmax(cman);

        setup_map_lat(cman, *cell_map1, *cell_map2);
	}
}
