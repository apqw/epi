#include <stdio.h>
#include <stdlib.h>

#include "global_cuda.h"
#include "cell_init.h"
#include "cell_connection_cuda.h"
#include "cell_interaction_cuda.h"
#include "map.h"
#include <iomanip>
#include <thread>

#include "ext_stim.h"
#include "cuda_host_util.h"
void output(DeviceData*d,int namei){
	static cell_pos_set pos[MAX_CELL_NUM];
	using namespace std;
	printf("memcpy start\n");
	cudaMemcpy(pos, d->c_pos_d[d->current], sizeof(cell_pos_set)*MAX_CELL_NUM, cudaMemcpyDeviceToHost);
	printf("memcpy end\n");
	//need last process end
	std::thread t1([namei,d](cell_pos_set* _pos){
		std::ofstream ofs(std::to_string(namei));
		for (int i = 0; i < d->ncell; i++){
			/*
			wfile << c->index << " "
			<< c->state << " "
			<< fixed << setprecision(15) << c->radius << " "
			<< fixed << setprecision(15) << c->ageb << " "
			<< fixed << setprecision(15) << c->agek << " "
			<< fixed << setprecision(15) << c->ca2p() << " "
			<< fixed << setprecision(15) << c->x() << " "
			<< fixed << setprecision(15) << c->y() << " "
			<< fixed << setprecision(15) << c->z() << " "
			<< fixed << setprecision(15) << c->ca2p_avg << " "
			<< c->rest_div_times << " "
			<< fixed << setprecision(15) << c->ex_fat << " "
			<< fixed << setprecision(15) << c->in_fat << " "
			<< (c->is_touch ? 1 : 0) << " "
			<< fixed << setprecision(15) << c->spr_nat_len << " "
			<< (int)(c->pair == nullptr ? (int)-1 : (int)(c->pair->index)) << " "
			<< c->fix_origin << std::endl;
			*/
			cell_pos_set cme = pos[i];
			unsigned int state = *(unsigned int*)&cme.w;
			//test
			ofs << i << " "
				<< state << " "
				<< fixed << setprecision(7) << (state == MEMB ? 1.0 : 1.4) << " "
				<< fixed << setprecision(7) << 0 << " "
				<< fixed << setprecision(7) << 0 << " "
				<< fixed << setprecision(7) << 0.122 << " "
				<< fixed << setprecision(7) << cme.x << " "
				<< fixed << setprecision(7) << cme.y << " "
				<< fixed << setprecision(7) << cme.z << " "
				<< fixed << setprecision(7) << 0.122 << " "
				<< 0 << " "
				<< fixed << setprecision(7) << 0 << " "
				<< fixed << setprecision(7) << 0 << " "
				<< 0 << " "
				<< fixed << setprecision(7) << 0 << " "
				<< -1 << " "
				<< 0 << std::endl;
		}
		return;
	},pos);
	t1.detach();
	
}


int main(int argc,char** argv)
{
	cudaFree(0);
	cudaDeviceSetLimit(cudaLimitMallocHeapSize, 1073741824);
	DeviceData a;
	init_with_file(&a, argv[1]);
	printf("%d\n", sizeof(connected_index_set)*50000/(1024*1024));
	
	size_t sz;
	cudaDeviceGetLimit(&sz, cudaLimitMallocHeapSize);
	printf("%lu malloc limit\n", sz);
	int* cmap1; float* cmap2,*ext_stim1,*ext_stim2;
	cudaMalloc((void**)&cmap1, sizeof(int)*(NX + 1)*(NY + 1)*(NZ + 1));
	cudaMalloc((void**)&cmap2, sizeof(float)*(NX + 1)*(NY + 1)*(NZ + 1));
	cudaMalloc((void**)&ext_stim1, sizeof(float)*(NX + 1)*(NY + 1)*(NZ + 1));
	cudaMalloc((void**)&ext_stim2, sizeof(float)*(NX + 1)*(NY + 1)*(NZ + 1));
	float* ext_stim_set[2] = { ext_stim1, ext_stim2 };
	cudaError_t ee;
	map_calc_init();
	connect_cell(&a);
	cell_pos_set zzmax_tmp;

	//thrust::device_ptr<cell_pos_set> pos_dptr[2] = { thrust::device_pointer_cast(a.c_pos_d[0]), thrust::device_pointer_cast(a.c_pos_d[1]) };

	for (int i = 0; i < 4000000; i++){
		/*
		if ((ee = cudaGetLastError()) != 0){
			printf("%d error\n", ee);
			system("pause");
			exit(1);
		}
		*/
		//system("pause");
		interact(&a);
		a.current = 1 - a.current;
		if (i % 100==0)printf("loop:%d done.\n", i);
		if (i % 10 == 0){
			
			connect_cell(&a);
		}
		
		setup_map(&a, cmap1, cmap2);
		
		//float aa = calc_zzmax_2(&a);
		//if (i % 100 == 0)printf("%lf test\n", aa);
		calc_ext_stim(&a, ext_stim_set[a.current], ext_stim_set[1 - a.current], cmap1, cmap2);
		if (i % 1000==0)output(&a,i/1000);
	}
	printf("done??\n");
	cudaFree(cmap1);
	cudaFree(cmap2);
	//init_grid(a.ncell,a.c_state_d,a.c_pos_d,a.c_connected_index_d)
	//a.~DeviceData();
	system("pause");
}