#include "global_cuda.h"
#include <string>
#include <fstream>
#include <cassert>
#include "cell_init.h"

template<unsigned int d>
int get_adj_memb_idx(int my_memb_idx){
	return 0;
}//dummy

template<>
int get_adj_memb_idx<(unsigned int)U>(int my_memb_idx){
	int kk= my_memb_idx / NMX;
	return ((kk + 1) % NMY)*NMX + my_memb_idx%NMX;
}

template<>
int get_adj_memb_idx<(unsigned int)L>(int my_memb_idx){
	int jj = my_memb_idx%NMX;
	return ((int)(my_memb_idx / NMX))*NMX + (jj - 1+NMX) % NMX;
}

template<>
int get_adj_memb_idx<(unsigned int)B>(int my_memb_idx){
	int kk = my_memb_idx / NMX;
	return ((kk - 1+NMY) % NMY)*NMX + my_memb_idx%NMX;
}

template<>
int get_adj_memb_idx<(unsigned int)R>(int my_memb_idx){
	int jj = my_memb_idx%NMX;
	return ((int)(my_memb_idx / NMX))*NMX + (jj + 1) % NMX;
}

void memb_init(int _nmemb, connected_index_set* cs){
	for (int j = 0; j<_nmemb; j++){
		const size_t jj = j%NMX;
		const size_t kk = j / NMX;
		cs[j].index[0] = get_adj_memb_idx<U>(j);
		cs[j].index[1] = get_adj_memb_idx<L>(j);
		cs[j].index[2] = get_adj_memb_idx<B>(j);
		cs[j].index[3] = get_adj_memb_idx<R>(j);
		//assert(cs[j].connect_index[0] >= 0 && cs[j].connect_index[1] >= 0 && cs[j].connect_index[2] >= 0 && cs[j].connect_index[3] >= 0);
	}
	//2-pass
	for (int j = 0; j < _nmemb; j++){

		const int memb_u = cs[j].index[0];
		const int memb_l = cs[j].index[1];
		const int memb_b = cs[j].index[2];
		const int memb_r = cs[j].index[3];

		//auto& memb_lu = cs[j].index[4];
		//auto& memb_bl = cs[j].index[5];
		//auto& memb_rb = cs[j].index[6];
		//auto& memb_ur = cs[j].index[7];

		cs[j].index[4] = cs[memb_l].index[0];
		cs[j].index[5] = cs[memb_b].index[1];
		cs[j].index[6] = cs[memb_r].index[2];
		cs[j].index[7] = cs[memb_u].index[3];

		cs[j].connected_num = 8;
		//assert(cs[j].connect_index[4] >= 0 && cs[j].connect_index[5] >= 0 && cs[j].connect_index[6] >= 0 && cs[j].connect_index[7] >= 0);
	}

}

void init_with_file(DeviceData* d, const char* path)
{
	/*
	ファイル読み込み試行
	*/
	std::ifstream dstrm;
	dstrm.open(path, std::ios::in);

	if (!dstrm) {
		assert(dstrm);
		throw "load failed";
	}

	using namespace cont;


	/*
	読み込み用一時変数
	*/
	CELL_STATE* c_state_h=new CELL_STATE[MAX_CELL_NUM];
	cell_pos_set* c_pos_h=new cell_pos_set[MAX_CELL_NUM];
	int* c_fix_origin_h = new int[MAX_CELL_NUM];
	connected_index_set* c_connected_index_h=new connected_index_set[MAX_CELL_NUM];
	int* c_pair_index_h = new int[MAX_CELL_NUM];
	float* c_ca2p_h = new float[MAX_CELL_NUM];
	float* c_ca2p_avg_h = new float[MAX_CELL_NUM];
	float* c_ex_inert_h = new float[MAX_CELL_NUM];
	float* c_agek_h = new float[MAX_CELL_NUM];
	float* c_ageb_h = new float[MAX_CELL_NUM];
	float* c_ex_fat_h = new float[MAX_CELL_NUM];
	float* c_in_fat_h = new float[MAX_CELL_NUM];
	float* c_spr_nat_len_h = new float[MAX_CELL_NUM];
	float* c_radius_h = new float[MAX_CELL_NUM];
	int* c_rest_div_times_h = new int[MAX_CELL_NUM];
	//int* c_dermis_index_h = new int[MAX_CELL_NUM];
	//CELL_STATE state = UNUSED;
	std::string line;
	int id_count = 0;

	unsigned int phase = 0;
	int nmemb = 0;
	int nder = 0;

	while (std::getline(dstrm, line)) {

		sscanf(line.c_str(), "%*d %d %f %f %f %f %f %f %f %f %d %f %f %*d %f %d %d",
			&c_state_h[id_count],
			&c_radius_h[id_count],
			&c_ageb_h[id_count],
			&c_agek_h[id_count],
			&c_ca2p_h[id_count],
			&c_pos_h[id_count].x,
			&c_pos_h[id_count].y,
			&c_pos_h[id_count].z,
			&c_ca2p_avg_h[id_count],
			&c_rest_div_times_h[id_count],
			&c_ex_fat_h[id_count],
			&c_in_fat_h[id_count],
			&c_spr_nat_len_h[id_count],
			&c_pair_index_h[id_count],
			&c_fix_origin_h[id_count]);

		/*
		BLANKにたどり着いたら終了
		*/
		CELL_STATE state = c_state_h[id_count];
		float rad = c_radius_h[id_count];
		if (state == BLANK)break;

		/*
		validation
		*/
		if (SYSTEM == BASAL && (state == ALIVE || state == DEAD || state == AIR)) {
			printf(" input date must not contain ALIVE or DEAD in case of BASAL\n");
			exit(1);
		}
		if (state == DER && rad != R_der) {
			printf("radii of DER not consistent with param.h %lf %lf\n",rad,R_der);
			exit(1);
		}
		if (state == MEMB && rad != R_memb) {
			printf("radii of DER not consistent with param.h\n");
			exit(1);
		}
		if (phase == 0 && state != MEMB) {
			assert(state == DER);
			phase++;
		}

		if (phase == 1 && state != DER) {
			//assert(state == DER);
			phase++;
		}
		if (phase > 0 && state == MEMB) {
			printf("non phase0 memb\n");
			exit(1);
		}
		if (phase > 1 && state == DER) {
			printf("non phase1 der\n");
			exit(1);
		}

		//if (state == FIX)printf("FIX\n");

		if (state == MEMB)nmemb++;
		if (state == DER)nder++;
		*(unsigned int*)(&c_pos_h[id_count].w) = c_state_h[id_count];
		printf("Phase %d  Cell loaded:%d\n", phase, id_count++);

	}

	assert(nmemb == NMX*NMY);
	d->ncell = id_count;
	d->nmemb = nmemb;
	d->nder = nder;
	d->nother = id_count - (nmemb + nder);

	memb_init(nmemb, c_connected_index_h);
	/*
		&c_state_h[id_count],
			&c_radius_h[id_count],
			&c_ageb_h[id_count],
			&c_agek_h[id_count],
			&c_ca2p_h[id_count],
			&c_pos_h[id_count].x,
			&c_pos_h[id_count].y,
			&c_pos_h[id_count].z,
			&c_ca2p_avg_h[id_count],
			&c_rest_div_times_h[id_count],
			&c_ex_fat_h[id_count],
			&c_in_fat_h[id_count],
			&c_spr_nat_len_h[id_count],
			&c_pair_index_h[id_count],
			&c_fix_origin_h[id_count]);
	*/
	cudaMemcpy(d->c_state_d, c_state_h, sizeof(CELL_STATE)*MAX_CELL_NUM, cudaMemcpyHostToDevice);
	//cudaMemcpy(d->c_radius_d, c_radius_h, sizeof(float)*MAX_CELL_NUM, cudaMemcpyHostToDevice);
	cudaMemcpy(d->c_ageb_d, c_ageb_h, sizeof(float)*MAX_CELL_NUM, cudaMemcpyHostToDevice);
	cudaMemcpy(d->c_agek_d, c_agek_h, sizeof(float)*MAX_CELL_NUM, cudaMemcpyHostToDevice);
	cudaMemcpy(d->c_ca2p_d[0], c_ca2p_h, sizeof(float)*MAX_CELL_NUM, cudaMemcpyHostToDevice);
	cudaMemcpy(d->c_ca2p_d[1], c_ca2p_h, sizeof(float)*MAX_CELL_NUM, cudaMemcpyHostToDevice);
	cudaMemcpy(d->c_pos_d[0], c_pos_h, sizeof(cell_pos_set)*MAX_CELL_NUM, cudaMemcpyHostToDevice);
	cudaMemcpy(d->c_pos_d[1], c_pos_h, sizeof(cell_pos_set)*MAX_CELL_NUM, cudaMemcpyHostToDevice);
	cudaMemcpy(d->c_ca2p_avg_d, c_ca2p_avg_h, sizeof(float)*MAX_CELL_NUM, cudaMemcpyHostToDevice);
	cudaMemcpy(d->c_rest_div_times_d, c_rest_div_times_h, sizeof(int)*MAX_CELL_NUM, cudaMemcpyHostToDevice);
	cudaMemcpy(d->c_ex_fat_d, c_ex_fat_h, sizeof(float)*MAX_CELL_NUM, cudaMemcpyHostToDevice);
	cudaMemcpy(d->c_in_fat_d, c_in_fat_h, sizeof(float)*MAX_CELL_NUM, cudaMemcpyHostToDevice);
	cudaMemcpy(d->c_spr_nat_len_d, c_spr_nat_len_h, sizeof(float)*MAX_CELL_NUM, cudaMemcpyHostToDevice);
	cudaMemcpy(d->c_pair_index_d, c_pair_index_h, sizeof(int)*MAX_CELL_NUM, cudaMemcpyHostToDevice);
	cudaMemcpy(d->c_fix_origin_d, c_fix_origin_h, sizeof(int)*MAX_CELL_NUM, cudaMemcpyHostToDevice);
	cudaMemcpy(d->c_connected_index_d, c_connected_index_h, sizeof(connected_index_set)*MAX_CELL_NUM, cudaMemcpyHostToDevice);
}