#include "proc.h"
#include "swapdata.h"
#include "Field.h"
#include "cell.h"
#include "cellmanager.h"
#include "define.h"


void proc(std::string init_data_path,std::string param_path,std::string init_uvp_data,std::string init_ATP_data,std::string init_ext_stim_data)
{
	using namespace cont;
	init_param_loader();
	load_param(param_path);
	Cell::cells.clear();
	Cell::load_from_file(init_data_path);
	
	auto ATP=new SwapData<Field<NX, NY, NZ>>();
	auto ext_stim= new SwapData<Field<NX, NY, NZ>>();
	printf("current cell num:%d\n", Cell::cells.size());

}
