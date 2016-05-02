#pragma once
#include <memory>
#include <fstream>
#include <thread>
#include <array>
#include "cell.h"
#include "const.h"
#include "static_if.hpp"
#include <atomic>
#define cdef static constexpr

class Field {
private:
	/*
		threads
	*/
	std::array<std::thread, THREAD_NUM> threads;
	CELL_STATE th_cst;
	double th_cx, th_cy, th_cz, th_crad;
	CellPtr th_RO_cell;
	std::vector<Arr2D<CellPtr>> th_out_conn_cell;
	std::array<int, THREAD_NUM> th_out_conn_count = { 0 };
	std::array<bool, THREAD_NUM> th_out_touch = { false };
	std::array<std::atomic<bool>, THREAD_NUM> th_need_calc;

	/*
		Cell management
	*/
	Arr1D<CellID> unused_id;
	int next_id_idx;
	void init_unused_id(Arr1D<CellID> &id_arr, int &next_id_var);
	int generate_new_id();
	void recycle_unused_id();
	void clean_up();
	void born_update();
	void state_update();
	void alloc_id_ptr(const Arr1D<CellPtr>& carr);
	void realloc_id_ptr_all();
	void update_all_cells();
	void copy_cell_data_to_old();
	void update_data();
	CellPtr add_new_cell_direct(const CELL_STATE& cstate, const CellData &c_data);

	/*
		Cell geom
	*/
	Arr1D<int> memb_u;
	Arr1D<int> memb_b;
	Arr1D<int> memb_l;
	Arr1D<int> memb_r;
	Arr1D<int> memb_ll;
	Arr1D<int> memb_bb;
	bool memb_indices_are_set;
	Arr3D<Arr2D<CellPtr>> cell_grid;
	Arr3D<int> cell_grid_count;
	void set_memb_indices();
	void cell_periodic_pos_fix();
	void update_cell_grid();
	void connect_cell_grid_ranged(int th_num, std::atomic<bool>* need_calc, CELL_STATE * cst, double * cx, double * cy, double * cz, double * crad, const CellPtr * read_only_cell, Arr2D<CellPtr>* out_conn_cell, int * out_conn_count, bool * out_touch);
	void connect_cell();
	void set_cell_lattice();
	double get_z_max();

	/*
		Cell dynamics
	*/
	void pair_interact();
	void pair_div();
	void memb_bending();
	void state_renew();
	void cell_dyn();

	/*
		Field mapping
	*/
	void setup_map();


	/*
		Ca dynamics
	*/
	class map_setter;
	class ca_dead_proc;
	class ca_common_cell_proc;
	class ca_common_cell_proc_conn;
	class ca_loop_init;
	class ca_result_setter;
	class ca_dyn_init_current;
	void ca_dyn();
	void calc_ext_stim();

	template<CELL_STATE my_state, CELL_STATE proc_state, bool use, class fobj>
	void __cell_interac_apply(CellPtr& me) {
		//STATIC_IF part will adaptively be deleted by compiler optimization
		auto& conn = me->old_data.connected_cell;
		STATIC_IF(use && my_state != proc_state) {
			for (auto &c : conn[proc_state]) {

				cell_intr_tmpl<fobj>(me, c);
			}
		}STATIC_ELSEIF(use && my_state == proc_state) {
			for (auto &c : conn[proc_state]) {
				if (me > c) {
					cell_intr_tmpl<fobj>(me, c);
				}
			}
		}STATIC_ENDIF;
	}

	template<CELL_STATE my_state,
		bool use_memb = false, class membf = null_coef,
		bool use_der = false, class derf = null_coef,
		bool use_alive = false, class alivef = null_coef,
		bool use_fix = false, class fixf = null_coef,
		bool use_musume = false, class musumef = null_coef,
		bool use_air = false, class airf = null_coef,
		bool use_dead = false, class deadf = null_coef
	>
		void cell_interact() {
		for (auto &me : a_cell[my_state]) {
			STATIC_IF(my_state == MEMB || my_state == DER) {
				if (me->old_data.z < me->old_data.radius) {
					interac_wall(me);
				}
			}STATIC_ENDIF;

			/*
			This part should be moved.
			*/
			STATIC_IF(my_state == FIX || my_state == MUSUME) {
				me->old_data.dermis = find_dermis(me);
				//me->old_data.spring_force = get_spring_force(me);
			}STATIC_ENDIF;

			__cell_interac_apply<my_state, MEMB, use_memb, membf>(me);
			__cell_interac_apply<my_state, DER, use_der, derf>(me);
			__cell_interac_apply<my_state, ALIVE, use_alive, alivef>(me);
			__cell_interac_apply<my_state, FIX, use_fix, fixf>(me);
			__cell_interac_apply<my_state, MUSUME, use_musume, musumef>(me);
			__cell_interac_apply<my_state, AIR, use_air, airf>(me);
			__cell_interac_apply<my_state, DEAD, use_dead, deadf>(me);

		}
	}

public:

	/*
		a_cell		:ç◊ñEëSÇƒ
		ca2p		:Ca2+îZìx
		ext_stim	:ç◊ñEäOéhåÉï®éøîZìx
		cell_by_id	:cell_by_id[id]Ç…ÇªÇÃidÇÃç◊ñEÇ÷ÇÃÉ|ÉCÉìÉ^Ç™Ç†ÇÈ
	*/

	Arr2D<CellPtr> a_cell;
	Arr3D<double> ATP; Arr3D<double> old_ATP;
	Arr3D<double> ext_stim; Arr3D<double> old_ext_stim;
	Arr3D<int> cornif_map; Arr3D<double> age_map;
	Arr3D<double> air_stim;
	Arr1D<CellPtr> cell_by_id;
	Arr1D<CellPtr> paired_cell_lower_id;
	Arr2D<CellPtr> pending_born_cell;
	Arr3D<CellPtr*> cell_map;
	Arr3D<double> cell_map2;

	int dead_count; //not disa
	double zzmax;
	int c_step = 0;
	bool is_sc_forced = false;
	int num_sc = 0;

	Field()
	{
		dead_count = 0;
		a_cell = InitArr2D<CellPtr>(STATE_NUM, 0,nullptr);
		ATP = InitArr3D<double>(NX + 1, NY + 1, NZ + 1,a0);
		old_ATP = InitArr3D<double>(NX + 1, NY + 1, NZ + 1, a0);
		ext_stim = InitArr3D<double>(NX + 1, NY + 1, NZ + 1, B0);
		old_ext_stim = InitArr3D<double>(NX + 1, NY + 1, NZ + 1, B0);
		air_stim = InitArr3D<double>(NX + 1, NY + 1, NZ + 1, 0);
		//ca2p_inert = Arr2D<double>(MAX_CELL_NUM + 1, Arr1D<double>(MAX_CELL_NUM + 1, 0));
		cell_by_id = Arr1D<CellPtr>(MAX_CELL_NUM + 1, nullptr);
		paired_cell_lower_id = Arr1D<CellPtr>(0);
		pending_born_cell = InitArr2D<CellPtr>(STATE_NUM, 0, nullptr);

		cell_grid = InitArr3D<Arr2D<CellPtr>>(ANX, ANY, ANZ, InitArr2D<CellPtr>(STATE_NUM, 0, nullptr));
		cell_grid_count = InitArr3D<int>(ANX, ANY, ANZ, 0);//Arr3D<int>(ANX, Arr2D<int>(ANY, Arr1D<int>(ANZ, 0)));
		th_out_conn_cell = InitArr3D<CellPtr>(THREAD_NUM, STATE_NUM, 0, nullptr);
		//do not change
		cell_map = InitArr3D<CellPtr*>(NX + 1, NY + 1, NZ + 1, nullptr);
		cell_map2 = InitArr3D<double>(NX + 1, NY + 1, NZ + 1, 0);
		age_map = InitArr3D<double>(NX + 1, NY + 1, NZ + 1, 0);
		cornif_map = InitArr3D<int>(NX + 1, NY + 1, NZ + 1, 0);

		init_unused_id(unused_id, next_id_idx);
	}

	void add_future_cell(const CELL_STATE& cstate, const CellData &c_data);
	void init_with_file(std::ifstream& data_stream);

	template<class CFunc, CELL_STATE first, CELL_STATE... cst>
	void proc_all_state_cell(Arr2D<CellPtr>& cell_set) {
		auto fobj = CFunc();
		int _end = cell_set[first].size();
		for (int i = 0; i < _end; i++) {
			fobj.operator() < first > (cell_set[first][i], this);
		}
		proc_all_state_cell<CFunc, cst...>(cell_set);
	}
	template<class CFunc>
	void proc_all_state_cell(Arr2D<CellPtr>& cell_set) {
	}

	template<class CFuncObj>
	void main_loop(CFuncObj callback_functor) {
		/*
			some init
		*/
		//is_sc_forced = true;
		genrand_init();
		/*
		for (int i = 0; i <= 4; i++) {
			th_need_calc[i] = false;
			 std::thread th([&,i] {
				connect_cell_grid_ranged(i,
					&(th_need_calc[i]),
					&th_cst, &th_cx, &th_cy, &th_cz, &th_crad, &th_RO_cell,
					&(th_out_conn_cell[i]), &(th_out_conn_count[i]), &(th_out_touch[i])
					);
			});
			 th.detach();
		}
		*/
		set_memb_indices();
		connect_cell();
		update_all_cells();

		//set_memb_indices();

		for (c_step = 0; c_step < NUM_ITR; c_step++) {
			std::cout << "LOOP:" << c_step << " MUSUME:" << a_cell[MUSUME].size() << " fix_age:" << a_cell[FIX][0]->data.ageb << std::endl;
			cell_dyn();
			update_all_cells();
			zzmax = get_z_max();
			if (SYSTEM == WHOLE) {
				set_cell_lattice();
				setup_map();
				calc_ext_stim();
				if (dead_count >= SW_THRESH || num_sc > 0) {
					printf("calc ca dyn\n");
					ca_dyn();
					if (num_sc > 0)num_sc--;
					dead_count = 0;
					printf("end\n");
				}
				update_all_cells();
			}
			//callback_functor(this);
		}
	}

};