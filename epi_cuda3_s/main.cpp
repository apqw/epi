#include "define.h"
#include "codetest.h"
#include "CData.h"
#include "CellDataMan.h"
#include <cstdlib>
#include "fsys.h"
#include <chrono>
void proc(const std::string& init_data_path, bool force_cornif,bool use_last = false, const std::string& init_uvp_data = "", const std::string& init_w_data = ""){
	CellDataMan cm;
	cm.init(init_data_path, use_last, init_uvp_data, init_w_data);

	make_dir(OUTPUT_DIR);
	int num_sc = 0;
	auto start = std::chrono::system_clock::now();
	cm.check_initialized();

}

int main(int argc,char** argv){
	codetest();
	proc(argv[1],true);
	pause();
}

