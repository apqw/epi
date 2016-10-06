#include "calc.h"
#include "../global.h"
#include "../define.h"
#include "cpu2/cell_connection.h"
#include "cpu2/CellManager.h"
#include "cpu2/cell_interaction.h"
#include "../fsys.h"
static void output_cell_data(CellManager&cman,int index){
    std::cout << "output:" << index << std::endl;
	cman.output(pm->outputdir+"/"+std::to_string(index));
}

static void initialize(CellManager&cman){
    cman._memb_init();
	connect_cell(cman);
	cman.init_value();
    make_dir(pm->outputdir.c_str());
}



static void cell_dynamics(CellManager&cman) {
    cell_interaction(cman);
    cman.pos_periodic_fix();
}
void calc_with_cpu(CellManager&cman){
    initialize(cman);

	for(size_t i=0;i<pm->NUM_ITR;i++){
		if(i%pm->CUT==0){
			output_cell_data(cman, static_cast<int>(i/pm->CUT));
		}
        cell_dynamics(cman);
        cman.pos_update();
	}
}
