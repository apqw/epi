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
#include "fsys.h"
#include <tbb/task_group.h>
#include <string>
#include <sstream>
#include <iostream>
#include <chrono>

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

void output_cell_data(CellManager& cman,size_t index){
    std::cout<<"logging..."<<std::endl;
    std::ostringstream _fname;
    _fname<<OUTPUTDIR<<"/"<<index;
    cman.output(_fname.str());
    std::cout<<"done."<<std::endl;
}

void output_uvp(CellManager& cman){
    FILE *fuvp;
    fuvp=std::fopen(last_data_uvp_name,"w");
    cman.all_foreach([fuvp](Cell* c){
        std::fprintf(fuvp,"%d %1.14e %1.14e %1.14e\n",(int)(c->get_index()),(double)(c->ca2p()),(double)(c->ex_inert),(double)(c->IP3()));
    });
    std::fclose(fuvp);
}

//need clean up
void output_w(CellManager& cman){
    FILE* fw_alt;
    fw_alt=std::fopen(last_data_w_name,"w");
    cman.all_foreach([fw_alt](Cell* c){
    std::fprintf(fw_alt,"%d %d\n",(int)(c->get_index()),(int)(c->gj.size()));
    for(auto& gjv:c->gj){
    std::fprintf(fw_alt,"%d %1.14e\n",(int)(gjv.first->get_index()),(double)(gjv.second()));
    }
    });
    std::fclose(fw_alt);
}

void output_field_data(CellManager& cman,const FArr3D<double>& ATP_first,const FArr3D<double>& ext_stim_first){
cman.output(last_data_cell_name);
output_uvp(cman);

output_w(cman);

ATP_first.output_binary(last_data_a_name);
ext_stim_first.output_binary(last_data_B_name);

}

void proc(const std::string& init_data_path,bool use_last,const std::string& init_uvp_data,const std::string& init_w_data,const std::string& init_ATP_data,
          const std::string& init_ext_stim_data)
{
	using namespace cont;
	auto cellset = std::make_unique<CellManager>();
    auto ATP=std::make_unique<SwapData<FArr3D<double>>>();
    auto ext_stim= std::make_unique<SwapData<FArr3D<double>>>();
    auto cell_map1 = std::make_unique<FArr3D<const Cell*>>();
    auto cell_map2 = std::make_unique<FArr3D<uint_fast8_t>>();
    cman_init(*cellset, init_data_path,use_last,init_uvp_data,init_w_data);
    if(use_last){
        ATP->first().read_binary(init_ATP_data);
        ATP->second().read_binary(init_ATP_data);
        ext_stim->first().read_binary(init_ext_stim_data);
        ext_stim->second().read_binary(init_ext_stim_data);
    }
	init_precalc_lat();
	init_precalc_per();
    make_dir(OUTPUTDIR);
	bool flg_forced_sc = FORCE_CORNIF;
	int num_sc = 0;

	double zzmax = 0;
    printf("current cell num:%zd\n", cellset->size());
	auto& cman = *cellset;

auto start=std::chrono::system_clock::now();
    for (size_t i = 0; i < NUM_ITR; i++) {
        if(i%100==0){

            auto dur=std::chrono::system_clock::now()-start;
            printf("loop:%d elapsed[sec]:%lf\n", i,0.001*std::chrono::duration_cast<std::chrono::milliseconds>(dur).count());
        }
        if(i%CUT==0){
            output_cell_data(cman,i/CUT);
        }
        if(i%(CUT*10)==0&&i>0){
            printf("cleaning up...\n");
            cman.clean_up();
            printf("done.\n");
            printf("saving field data...\n");
            output_field_data(cman,ATP->first(),ext_stim->first());
            printf("done.\n");
        }

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
