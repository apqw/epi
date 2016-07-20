/*
 * main.cpp
 *
 *  Created on: 2016/07/07
 *      Author: yasu7890v
 */

#include "code_test.h"
#include "stdio.h"
#include "CellManager.h"
#include "cell_connection.h"
#include "cell_interaction.h"
#include "cell_state_renew.h"
#include "ext_stim.h"
#include "map.h"
#include <cstdlib>
#include <vector>
#include <fstream>
#include <string>
#include "utils.h"
#include "ca2p.h"
#include <iostream>
#include <iomanip>
#include "fsys.h"

void dbg_output(CellManager* cm,int idx){
	static std::vector<real4> pos(MAX_CELL_NUM);
	cudaMemcpy(&pos[0], cm->current_pos_host(), sizeof(real4)*MAX_CELL_NUM, cudaMemcpyDeviceToHost);
	std::ofstream ofs(("dbg" + int_to_string(idx)).c_str());
	for (int i = 0; i < cm->ncell_host; i++){
		ofs << i<<" "<<pos[i].x << " " << pos[i].y << " " << pos[i].z << std::endl;
	}

}

void output(const std::string &filename,CellManager* cm)
{
	std::ofstream wfile(filename);
	
	if (!wfile){
		std::cout << "Output file creation error. Filename:" << filename << std::endl;
		system("pause");
		exit(1);
	}
	

		using namespace std;
		static vector<CELL_STATE> state(MAX_CELL_NUM);
		static vector<real> ageb(MAX_CELL_NUM);
		static vector<real> agek(MAX_CELL_NUM);
		static vector<real> ca2p(MAX_CELL_NUM);
		static vector<CellPos> pos(MAX_CELL_NUM);
		static vector<real> ca2p_avg(MAX_CELL_NUM);
		static vector<int> rest_div_times(MAX_CELL_NUM);
		static vector<real> ex_fat(MAX_CELL_NUM);
		static vector<real> in_fat(MAX_CELL_NUM);
		//static vector<real> ageb(MAX_CELL_NUM);
		static vector<real> spr_nat_len(MAX_CELL_NUM);
		static vector<CellIndex> pair_idx(MAX_CELL_NUM);
		static vector<int> fix_origin(MAX_CELL_NUM);
		cudaMemcpy(&state[0], cm->state, sizeof(CELL_STATE)*MAX_CELL_NUM, cudaMemcpyDeviceToHost);
		cudaMemcpy(&ageb[0], cm->ageb, sizeof(real)*MAX_CELL_NUM, cudaMemcpyDeviceToHost);
		cudaMemcpy(&agek[0], cm->agek, sizeof(real)*MAX_CELL_NUM, cudaMemcpyDeviceToHost);
		cudaMemcpy(&ca2p[0], cm->ca2p[cm->current_phase_host], sizeof(real)*MAX_CELL_NUM, cudaMemcpyDeviceToHost);
		cudaMemcpy(&pos[0], cm->current_pos_host(), sizeof(CellPos)*MAX_CELL_NUM, cudaMemcpyDeviceToHost);
		cudaMemcpy(&ca2p_avg[0], cm->ca2p_avg, sizeof(real)*MAX_CELL_NUM, cudaMemcpyDeviceToHost);
		cudaMemcpy(&rest_div_times[0], cm->rest_div_times, sizeof(int)*MAX_CELL_NUM, cudaMemcpyDeviceToHost);
		cudaMemcpy(&ex_fat[0], cm->ex_fat, sizeof(real)*MAX_CELL_NUM, cudaMemcpyDeviceToHost);
		cudaMemcpy(&in_fat[0], cm->in_fat, sizeof(real)*MAX_CELL_NUM, cudaMemcpyDeviceToHost);
		cudaMemcpy(&spr_nat_len[0], cm->spr_nat_len, sizeof(real)*MAX_CELL_NUM, cudaMemcpyDeviceToHost);
		cudaMemcpy(&pair_idx[0], cm->pair_index, sizeof(CellIndex)*MAX_CELL_NUM, cudaMemcpyDeviceToHost);
		cudaMemcpy(&fix_origin[0], cm->fix_origin, sizeof(int)*MAX_CELL_NUM, cudaMemcpyDeviceToHost);
		for (int i = 0; i < cm->ncell_host; i++){
			wfile << i << " "
				<< state[i] << " "
				<< fixed << setprecision(15) << (state[i]==MEMB?R_memb:R_max) << " "
				<< fixed << setprecision(15) << ageb[i] << " "
				<< fixed << setprecision(15) << agek[i] << " "
				<< fixed << setprecision(15) << ca2p[i] << " "
				<< fixed << setprecision(15) << pos[i].x << " "
				<< fixed << setprecision(15) << pos[i].y << " "
				<< fixed << setprecision(15) << pos[i].z << " "
				<< fixed << setprecision(15) << ca2p_avg[i] << " "
				<< rest_div_times[i] << " "
				<< fixed << setprecision(15) << ex_fat[i] << " "
				<< fixed << setprecision(15) << in_fat[i] << " "
				<< 0 << " " //touch
				<< fixed << setprecision(15) << spr_nat_len[i] << " "
				<< (int)(pair_idx[i] <0 ? (int)-1 : (int)(pair_idx[i])) << " "
				<< fix_origin[i] << std::endl;
		}
	
}

int main(int argc,char** argv){
	cudaDeviceSetLimit(cudaLimitMallocHeapSize, (size_t)1073741824);
	codetest();
	cudaError_t ee;
	bool forced_cornif = true;
	int num_sc = 0;
	make_dir("output");
	CellManager cm;
	if ((ee = cudaGetLastError()) != 0){
		printf("error ff:%d\n", ee);
		system("pause");
		exit(1);
	}
	cm.init_with_file(argv[1]);
	printf("loop start\n");
	
	//dbg_output(&cm, -1);
	connect_cell(&cm);
	int* cmap1; float* cmap2; real *ext_stim1, *ext_stim2;
	cudaMalloc((void**)&cmap1, sizeof(int)*(NX + 1)*(NY + 1)*(NZ + 1));
	cudaMalloc((void**)&cmap2, sizeof(float)*(NX + 1)*(NY + 1)*(NZ + 1));
	cudaMalloc((void**)&ext_stim1, sizeof(real)*(NX + 1)*(NY + 1)*(NZ + 1));
	cudaMalloc((void**)&ext_stim2, sizeof(real)*(NX + 1)*(NY + 1)*(NZ + 1));
	
	real* ext_stim_set[2] = { ext_stim1, ext_stim2 };
	map_calc_init();

	for(int i=0;i<NUM_ITR;i++){
		if ((ee = cudaGetLastError()) != 0){
			printf("error:%d\n", ee);
			cudaFree(cmap1);
			cudaFree(cmap2);
			cudaFree(ext_stim1);
			cudaFree(ext_stim2);
			system("pause");
			exit(1);
		}
		
		cell_interact(&cm);
		if ((ee = cudaGetLastError()) != 0){
			printf("error2:%d\n", ee);
			cudaFree(cmap1);
			cudaFree(cmap2);
			cudaFree(ext_stim1);
			cudaFree(ext_stim2);
			system("pause");
			exit(1);
		}
		cm.switch_phase();
		if ((ee = cudaGetLastError()) != 0){
			printf("error3:%d\n", ee);
			cudaFree(cmap1);
			cudaFree(cmap2);
			cudaFree(ext_stim1);
			cudaFree(ext_stim2);
			system("pause");
			exit(1);
		}
		cell_pos_periodic_fix(&cm);
		if ((ee = cudaGetLastError()) != 0){
			printf("error4:%d\n", ee);
			cudaFree(cmap1);
			cudaFree(cmap2);
			cudaFree(ext_stim1);
			cudaFree(ext_stim2);
			system("pause");
			exit(1);
		}
		cell_state_renew(&cm);
		if ((ee = cudaGetLastError()) != 0){
			printf("error5:%d\n", ee);
			cudaFree(cmap1);
			cudaFree(cmap2);
			cudaFree(ext_stim1);
			cudaFree(ext_stim2);
			system("pause");
			exit(1);
		}
		if (i % 1 == 0 || cm.should_force_reconnect())connect_cell(&cm);
		if ((ee = cudaGetLastError()) != 0){
			printf("error6:%d\n", ee);
			cudaFree(cmap1);
			cudaFree(cmap2);
			cudaFree(ext_stim1);
			cudaFree(ext_stim2);
			system("pause");
			exit(1);
		}
		setup_map(&cm, cmap1, cmap2);
		if ((ee = cudaGetLastError()) != 0){
			printf("error7:%d\n", ee);
			cudaFree(cmap1);
			cudaFree(cmap2);
			cudaFree(ext_stim1);
			cudaFree(ext_stim2);
			system("pause");
			exit(1);
		}
		real zmax = get_cell_zmax(&cm);
		if ((ee = cudaGetLastError()) != 0){
			printf("error8:%d\n", ee);
			cudaFree(cmap1);
			cudaFree(cmap2);
			cudaFree(ext_stim1);
			cudaFree(ext_stim2);
			system("pause");
			exit(1);
		}
		//printf("%lf zmax\n", zmax);
	
		calc_ext_stim(&cm, ext_stim_set[cm.current_phase_host], ext_stim_set[1 - cm.current_phase_host], cmap1, cmap2,zmax);
		if ((ee = cudaGetLastError()) != 0){
			printf("error9:%d\n", ee);
			cudaFree(cmap1);
			cudaFree(cmap2);
			cudaFree(ext_stim1);
			cudaFree(ext_stim2);
			system("pause");
			exit(1);
		}
		//printf("zmax:%f\n",zmax);
		if (i*DT_Cell > T_TURNOVER&&forced_cornif){
			forced_cornif = false;
			printf("forced cornif\n");
			initialize_sc(&cm, zmax);
			num_sc = NUM_SC_INIT;
		}
		if (cm.sw_host >= SW_THRESH || num_sc > 0){
			printf("calc ca...\n");
			calc_ca2p(&cm, ext_stim_set[cm.current_phase_host], cmap1, cmap2, zmax);
			printf("end.\n");
			if(num_sc>0)num_sc--;
			cudaMemset(cm.sw, 0, sizeof(int));
		}

		//dbg_output(&cm, i);
		if(i%100==0)printf("tesuya %d\n",i);
		if (i % 1000 == 0&&i!=0){
			output("output/" + std::to_string((unsigned long long)(i / 1000)), &cm);
		}
	}
	cudaFree(cmap1);
	cudaFree(cmap2);
	cudaFree(ext_stim1);
	cudaFree(ext_stim2);
#ifdef _WIN32
	system("pause");
#endif
}
