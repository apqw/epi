#include "calc.h"
#include "../global.h"
#include "../define.h"
#include "cpu2/cell_connection.h"
#include "cpu2/CellManager.h"
static void output_cell_data(CellManager&cman,int index){
	cman.output(pm->outputdir+"/"+std::to_string(index));
}

static void initialize(CellManager&cman){

	connect_cell(cman);
	cman.init_value();
}

void calc_with_cpu(CellManager&cman){
	for(size_t i=0;i<pm->NUM_ITR;i++){
		if(i%pm->CUT==0){
			output_cell_data(cman, static_cast<int>(i/pm->CUT));
		}

	}
}
