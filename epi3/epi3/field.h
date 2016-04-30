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
	Arr1D<CellID> unused_id;
	int next_id_idx;
	void init_unused_id(Arr1D<CellID> &id_arr, int &next_id_var);
	int generate_new_id();
	void recycle_unused_id();
	void alloc_id_ptr(const Arr1D<CellPtr>& carr);
	void realloc_id_ptr_all();
	Arr3D<Arr2D<CellPtr>> cell_grid;
	Arr3D<int> cell_grid_count;

	Arr1D<int> memb_u;
	Arr1D<int> memb_b;
	Arr1D<int> memb_l;
	Arr1D<int> memb_r;
	Arr1D<int> memb_ll;
	Arr1D<int> memb_bb;
	bool memb_indices_are_set;
	CellPtr add_new_cell_direct(const CELL_STATE& cstate, const CellData &c_data);
	std::array<std::thread, THREAD_NUM> threads;
	CELL_STATE th_cst;
	double th_cx, th_cy, th_cz,th_crad;
	CellPtr th_RO_cell;
	std::vector<Arr2D<CellPtr>> th_out_conn_cell;
	std::array<int, THREAD_NUM> th_out_conn_count = { 0 };
	std::array<bool, THREAD_NUM> th_out_touch = { false };

	std::array<std::atomic<bool>, THREAD_NUM> th_need_calc;

public:

	/*
		a_cell		:ç◊ñEëSÇƒ
		ca2p		:Ca2+îZìx
		ext_stim	:ç◊ñEäOéhåÉï®éøîZìx
		cell_by_id	:cell_by_id[id]Ç…ÇªÇÃidÇÃç◊ñEÇ÷ÇÃÉ|ÉCÉìÉ^Ç™Ç†ÇÈ
	*/

	Arr2D<CellPtr> a_cell;
	Arr3D<double> ca2p;
	Arr3D<double> ext_stim;
	Arr1D<CellPtr> cell_by_id;
	Arr1D<CellPtr> paired_cell_lower_id;
	Arr2D<CellPtr> pending_born_cell;
	
	int dead_count; //not disa
	double zzmax;
	//ïsäàê´å¯â 
	//ID*ID
	//Arr2D<double> ca2p_inert;

	/*
		unused_id	:åªç›ñ¢äÑìñÇÃIDÇäiî[(ÉTÉCÉYÇÕMAX_CELL_NUM)
		next_id		:éüäÑÇËìñÇƒÇÈIDÇéwÇ∑unused_idè„ÇÃà íu
		IDÇÕ1Å`MAX_CELL_NUMÇ∆Ç∑ÇÈÅBID:0ÇÕinvalid
		unused_idÇÃèâä˙èÛë‘ÇÕunused_id[i]=i+1;
	*/

	
	Field()
	{
		a_cell = Arr2D<CellPtr>(STATE_NUM, Arr1D<CellPtr>(0,nullptr));
		ca2p = Arr3D<double>(NX, Arr2D<double>(NY, Arr1D<double>(NZ, 0)));
		ext_stim = Arr3D<double>(NX, Arr2D<double>(NY, Arr1D<double>(NZ, 0)));
		//ca2p_inert = Arr2D<double>(MAX_CELL_NUM + 1, Arr1D<double>(MAX_CELL_NUM + 1, 0));
		cell_by_id = Arr1D<CellPtr>(MAX_CELL_NUM+1, nullptr);
		paired_cell_lower_id = Arr1D<CellPtr>(0);
		pending_born_cell = Arr2D<CellPtr>(STATE_NUM, Arr1D<CellPtr>(0));

		cell_grid =
			Arr3D<Arr2D<CellPtr>>(ANX,
				Arr2D<Arr2D<CellPtr>>(ANY,
					Arr1D<Arr2D<CellPtr>>(ANZ,
						Arr2D<CellPtr>(STATE_NUM, Arr1D<CellPtr>(0)))));
		cell_grid_count = Arr3D<int>(ANX, Arr2D<int>(ANY, Arr1D<int>(ANZ, 0)));
			th_out_conn_cell = std::vector<Arr2D<CellPtr>>(THREAD_NUM,Arr2D<CellPtr>(STATE_NUM, Arr1D<CellPtr>(0, nullptr)));

	
		init_unused_id(unused_id, next_id_idx);
	}

	void add_future_cell(const CELL_STATE& cstate,const CellData &c_data);
	
	void clean_up();
	void born_update();
	void state_update();
	void update_all_cells();
	void copy_cell_data_to_old();
	void state_renew();
	void cell_dyn();
	void cell_periodic_pos_fix();
	void update_cell_grid();
	void connect_cell_grid_ranged(int th_num, std::atomic<bool>* need_calc, CELL_STATE * cst, double * cx, double * cy, double * cz, double * crad, const CellPtr * read_only_cell, Arr2D<CellPtr>* out_conn_cell, int * out_conn_count, bool * out_touch);
	void connect_cell();
	void update_data();
	void set_memb_indices();
	void pair_interact();
	void pair_div();
	void memb_bending();
	void calc_energy();
	void init_with_file(std::ifstream& data_stream);

	void ca_dyn();
	double get_z_max();


	template<CELL_STATE my_state,
		bool use_memb	= false,	class membf		= null_coef,
		bool use_der	= false,	class derf		= null_coef,
		bool use_alive	= false,	class alivef	= null_coef,
		bool use_fix	= false,	class fixf		= null_coef,
		bool use_musume = false,	class musumef	= null_coef,
		bool use_air	= false,	class airf		= null_coef,
		bool use_dead	= false,	class deadf		= null_coef
	>
	void cell_interact() {
		//memb interact

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


			auto& conn = me->old_data.connected_cell; //use ref
			//STATIC_IF part will adaptively be deleted by compiler optimization
			STATIC_IF(use_memb && my_state!=MEMB) {
				for (auto &cmemb : conn[MEMB]) {
	
						cell_intr_tmpl<membf>(me, cmemb);
				}
			}STATIC_ELSEIF(use_memb && my_state == MEMB) {
				for (auto &cmemb : conn[MEMB]) {
					if (me > cmemb) {
						cell_intr_tmpl<membf>(me, cmemb);
					}
				}
			}STATIC_ENDIF;


			STATIC_IF(use_der && my_state != DER) {
				for (auto &cder : conn[DER]) {
						cell_intr_tmpl<derf>(me, cder);
				}
			}STATIC_ELSEIF(use_der && my_state == DER) {
				for (auto &cder : conn[DER]) {
					if (me > cder) {
						cell_intr_tmpl<derf>(me, cder);
					}
				}
			}STATIC_ENDIF;



			STATIC_IF(use_alive && my_state!=ALIVE) {
				for (auto &calive : conn[ALIVE]) {
						cell_intr_tmpl<alivef>(me, calive);
				}
			}STATIC_ELSEIF(use_alive && my_state == ALIVE) {
				for (auto &calive : conn[ALIVE]) {
					if(me>calive)
					cell_intr_tmpl<alivef>(me, calive);
				}
			}STATIC_ENDIF;

			STATIC_IF(use_fix && my_state!=FIX) {
				
				for (auto &cfix : conn[FIX]) {
						cell_intr_tmpl<fixf>(me, cfix);
				}
			}STATIC_ELSEIF(use_fix && my_state == FIX) {
				for (auto &cfix : conn[FIX]) {
					if(me>cfix)
					cell_intr_tmpl<fixf>(me, cfix);
				}
			}STATIC_ENDIF;

			STATIC_IF(use_musume && my_state!=MUSUME) {
				for (auto &cmusume : conn[MUSUME]) {
						cell_intr_tmpl<musumef>(me, cmusume);
				}
			}STATIC_ELSEIF(use_musume && my_state == MUSUME) {
				for (auto &cmusume : conn[MUSUME]) {
					if(me>cmusume)cell_intr_tmpl<musumef>(me, cmusume);
				}
			}STATIC_ENDIF;

			STATIC_IF(use_air && my_state!=AIR) {
				for (auto &cair : conn[AIR]) {
						cell_intr_tmpl<airf>(me, cair);
				}
			}STATIC_ELSEIF(use_air && my_state == AIR) {
				for (auto &cair : conn[AIR]) {
					if(me>cair)cell_intr_tmpl<airf>(me, cair);
				}
			}STATIC_ENDIF;

			STATIC_IF(use_dead && my_state!=DEAD) {
				for (auto &cdead : conn[DEAD]) {
						cell_intr_tmpl<deadf>(me, cdead);
				}
			}STATIC_ELSEIF(use_dead && my_state == DEAD) {
				for (auto &cdead : conn[DEAD]) {
					if(me>cdead)cell_intr_tmpl<deadf>(me, cdead);
				}
			}STATIC_ENDIF;



		}
	}
	template<class CFunc, CELL_STATE first,CELL_STATE... cst>
	void proc_all_state_cell() {
		for (auto&c : a_cell[first]) {
			CFunc()(c);
		}
		proc_all_state_cell<CFunc, cst...>();
	}
	template<class CFunc>
	void proc_all_state_cell() {
	}


	template<class CFuncObj>
	void main_loop(CFuncObj callback_functor) {

		/*
			some init
			
		*/
		genrand_init();
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
		set_memb_indices();
		connect_cell();
		update_all_cells();

		//set_memb_indices();

		for (int i = 0; i < NUM_ITR; i++) {
			std::cout << "LOOP:" << i <<" MUSUME:"<< a_cell[MUSUME].size()<< std::endl;
			cell_dyn();
			update_all_cells();
			zzmax = get_z_max();
			//callback_functor(this);
		}
	}

};