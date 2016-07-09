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
#include <cstdlib>
#include <vector>
#include <fstream>
#include <string>

void dbg_output(CellManager* cm,int idx){
	static std::vector<float4> pos(MAX_CELL_NUM);
	cudaMemcpy(&pos[0], cm->current_pos_host(), sizeof(float4)*MAX_CELL_NUM, cudaMemcpyDeviceToHost);
	std::ofstream ofs("dbg" + std::to_string(idx));
	for (int i = 0; i < cm->ncell_host; i++){
		ofs << i<<" "<<pos[i].x << " " << pos[i].y << " " << pos[i].z << std::endl;
	}

}
int main(int argc,char** argv){
	cudaDeviceSetLimit(cudaLimitMallocHeapSize, 1073741824);
	codetest();

	CellManager cm;
	cm.init_with_file(argv[1]);
	printf("loop start\n");
	cudaError_t ee;
	//dbg_output(&cm, -1);
	connect_cell(&cm);
	for(int i=0;i<100000;i++){
		if ((ee = cudaGetLastError()) != 0){
			printf("error:%d\n", ee);
			system("pause");
			exit(1);
		}
		
		cell_interact(&cm);
		cm.switch_phase();
		cell_pos_periodic_fix(&cm);
		cell_state_renew(&cm);
		
		if (i % 10 == 0||cm.should_force_reconnect())connect_cell(&cm);
		//dbg_output(&cm, i);
		if(i%100==0)printf("tesuya %d\n",i);
	}
	system("pause");
}
