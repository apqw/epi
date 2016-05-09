#pragma once
#include <list>
#include <array>
#include <vector>
#include <map>
#include <cmath>
#include "define.h"
#include "Cell.h"
#include "util_func.h"



	using CellSpace = std::list<Cell*>;

	template<unsigned int DIV>
	class CellOctree {
		const double spatial_div_x = cont::LX / spow<int>(2, DIV);
		const double spatial_div_y = cont::LY / spow<int>(2, DIV);
		const double spatial_div_z = cont::LZ / spow<int>(2, DIV);


		const int snum =( spow<int>(8, DIV+1) - 1) / 7;
		std::array<CellSpace, (spow<int>(8, DIV + 1) - 1) / 7> slist;
		std::vector<Cell*> collide_list;
		std::vector<std::pair<Cell*, Cell*>> collide_candidate;

		void _collision(int space_index) {
			CellSpace& cur_sp = slist[space_index];
			for (auto& cptr : cur_sp) {
				for (auto& cptr2 : cur_sp) {
					if (cptr >= cptr2) {
						collide_candidate.push_back(std::make_pair(cptr, cptr2));
					}
				}
				for (auto& stacked_cell : collide_list) {
					collide_candidate.push_back(std::make_pair(cptr, stacked_cell));
				}
			}
			
			int child_head_index = space_index * 8 + 1;

			if (child_head_index >= snum) { //->no child
				return;
			}

			for (auto& cptr : cur_sp) {
				collide_list.push_back(cptr);
			}

			for (int i = 0; i < 8; i++) {
				int child_index = space_index * 8 + 1 + i;
				_collision(child_index);
			}

			for (auto& _dummy : cur_sp) {
				collide_list.pop_back();
			}

		}
		void set_cellspace(CellMan& cman) {
			using namespace cont;
			cman.foreach([&](CellPtr& c, int i) {
				double lg_x = c->pos[0]() + c->radius()*LJ_THRESH / 2.0;
				double lg_y = c->pos[1]() + c->radius()*LJ_THRESH / 2.0;
				double lg_z = c->pos[2]() + c->radius()*LJ_THRESH / 2.0;
				lg_z = lg_z > 0 ? lg_z : 0;

				double le_x = c->pos[0]() - c->radius()*LJ_THRESH / 2.0;
				double le_y = c->pos[1]() - c->radius()*LJ_THRESH / 2.0;
				double le_z = c->pos[2]() - c->radius()*LJ_THRESH / 2.0;
				le_z = le_z > 0 ? le_z : 0;


				uint32_t lg_ix = (uint32_t)((0.5*LX - p_diff_sc_x(0.5*LX, lg_x)) / spatial_div_x);
				uint32_t lg_iy = (uint32_t)((0.5*LY - p_diff_sc_y(0.5*LY, lg_y)) / spatial_div_y);
				uint32_t lg_iz = (uint32_t)(lg_z / spatial_div_z);
				uint32_t lg_m_order = Get3DMortonOrder(lg_ix, lg_iy, lg_iz);

				uint32_t le_ix = (uint32_t)((0.5*LX - p_diff_sc_x(0.5*LX, le_x)) / spatial_div_x);
				uint32_t le_iy = (uint32_t)((0.5*LY - p_diff_sc_y(0.5*LY, le_y)) / spatial_div_y);
				uint32_t le_iz = (uint32_t)(le_z / spatial_div_z);

				uint32_t le_m_order = Get3DMortonOrder(le_ix, le_iy, le_iz);

				uint32_t m_flag = lg_m_order^le_m_order;//XOR

				uint32_t shared_level = 0;
				uint32_t idx_in_shared = 0;
				for (unsigned int i = DIV - 1; i >= 0; i--) {
					if (((((unsigned)0b111) << (i * 3))&m_flag) != 0) {
						idx_in_shared = lg_m_order >> ((i + 1) * 3);
						break;
					}
					shared_level++;
				}

				slist[(spow<int>(8, shared_level) - 1) / 7 + idx_in_shared].push_back(c.get());

			});
		}
	public:
		template<class ColTestLambda>
		void collide_candidate_check(CellMan& cman,ColTestLambda&& ct) {
			for (auto& cs : slist) {
				cs.clear();
			}
			assert(collide_list.size() == 0);
			//collide_list.clear();
			collide_candidate.clear();
			set_cellspace(cman);
			_collision(0);
			for (auto& cpair : collide_candidate) {
				ct(cpair);
			}
		}
	};