#pragma once
#include "define.h"
#include "component.h"
#include "Cell.h"
#include "particle_seacrh.h"
#include <istream>
#include <fstream>
class Field
{
private:

	const int _MAX_CELL_NUM = 30000;

	CellMan cells;
	double _ATP[cont::NX + 1][cont::NY + 1][cont::NZ + 1];
	double old_ATP[cont::NX+1][cont::NY+1][cont::NZ+1];

	double _ext_stim[cont::NX + 1][cont::NY + 1][cont::NZ + 1];
	double old_ext_stim[cont::NX + 1][cont::NY + 1][cont::NZ + 1];
	Cell* cell_map[cont::NX + 1][cont::NY + 1][cont::NZ + 1];
	int cell_map2[cont::NX + 1][cont::NY + 1][cont::NZ + 1];
	int air_stim_flg[cont::NX + 1][cont::NY + 1][cont::NZ + 1];
    double* cell_diffu_map[cont::NX + 1][cont::NY + 1][cont::NZ + 1];
	void interact_cell();
	void cell_state_renew();
	void cell_pos_periodic_fix();
	void connect_cells();
	void set_cell_lattice();
	double calc_zzmax();
	double zzmax=0;
    int per_x_next_idx[cont::NX],per_x_prev_idx[cont::NX],per_y_next_idx[cont::NY],per_y_prev_idx[cont::NY];
    int sw=0;
    int num_sc=0;
    bool flg_forced_sc=false;
public:
	Field(
        int __MAX_CELL_NUM=30000,
            bool _forced_sc=false
		):_MAX_CELL_NUM(__MAX_CELL_NUM),
		_ATP(),
		old_ATP(),
		_ext_stim(),
		old_ext_stim(),
		cell_map(),
		cell_map2(),
        air_stim_flg(),
        cell_diffu_map(),flg_forced_sc(_forced_sc)
		{
using namespace cont;
std::memset(_ATP, 0, sizeof(double)*(NX + 1)*(NY + 1)*(NZ + 1));
std::memset(old_ATP, 0, sizeof(double)*(NX + 1)*(NY + 1)*(NZ + 1));
std::memset(_ext_stim, 0, sizeof(double)*(NX + 1)*(NY + 1)*(NZ + 1));
std::memset(old_ext_stim, 0, sizeof(double)*(NX + 1)*(NY + 1)*(NZ + 1));
std::memset(cell_map, 0, sizeof(Cell*)*(NX + 1)*(NY + 1)*(NZ + 1));
std::memset(cell_map2, 0, sizeof(int)*(NX + 1)*(NY + 1)*(NZ + 1));
std::memset(air_stim_flg, 0, sizeof(int)*(NX + 1)*(NY + 1)*(NZ + 1));
std::memset(cell_diffu_map, 0, sizeof(double*)*(NX + 1)*(NY + 1)*(NZ + 1));
        for(int k=0;k<NY;k++){
            int prev_y = k - 1;
            int next_y = k + 1;
            if (prev_y < 0) {
                prev_y += NY;
            }
            if (next_y >= NY) {
                next_y -= NY;
            }
            per_y_prev_idx[k]=prev_y;
            per_y_next_idx[k]=next_y;
        }

        for(int k=0;k<NX;k++){
            int prev_x= k - 1;
            int next_x = k + 1;
            if (prev_x < 0) {
                prev_x += NX;
            }
            if (next_x >= NX) {
                next_x -= NX;
            }
            per_x_prev_idx[k]=prev_x;
            per_x_next_idx[k]=next_x;
        }
	}


	void init_with_file(std::ifstream& dstrm);
	void cell_dynamics();
	void main_loop();
	void setup_map();
	void calc_b();
	void b_update();
	void calc_ca();
	void ATP_update();
    void initialize_sc();
    void check_localization();
    void output_data(int i);

	
	//~Field();
};

