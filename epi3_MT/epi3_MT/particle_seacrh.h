#pragma once
#include <list>
#include <array>
#include <vector>
#include <map>
#include <cmath>
#include "define.h"
#include "Cell.h"
#include "util_func.h"



	using CellSpace = std::vector<Cell*>;

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
			printf("col call1\n");
			CellSpace& cur_sp = slist[space_index];
			printf("cur size:%d\n", cur_sp.size());
			for (auto& cptr : cur_sp) {
				for (auto& cptr2 : cur_sp) {
					if (cptr > cptr2) {
						collide_candidate.push_back(std::make_pair(cptr, cptr2));
					}
				}
				for (auto& stacked_cell : collide_list) {
					collide_candidate.push_back(std::make_pair(cptr, stacked_cell));
				}
			}
			printf("col call2\n");
			int child_head_index = space_index * 8 + 1;

			if (child_head_index >= snum) { //->no child
				return;
			}
			printf("col call3\n");
			for (auto& cptr : cur_sp) {
				collide_list.push_back(cptr);
			}
			printf("col call4\n");
			for (int i = 0; i < 8; i++) {
				int child_index = space_index * 8 + 1 + i;
				_collision(child_index);
			}
			printf("col call5\n");
			for (auto& _dummy : cur_sp) {
				collide_list.pop_back();
			}
			printf("col call6\n");
		}
		void set_cellspace(CellMan& cman) {
			using namespace cont;
			cman.foreach([&](CellPtr& c, int cidx) {
				
				double lg_x = c->pos[0]() + c->radius()*LJ_THRESH;
				if (lg_x > LX) {
					lg_x -= LX;
				}
				else if (lg_x < 0) {
					lg_x += LX;
				}
				double lg_y = c->pos[1]() + c->radius()*LJ_THRESH;
				if (lg_y > LY) {
					lg_y -= LY;
				}
				else if (lg_y < 0) {
					lg_y += LY;
				}
				double lg_z = c->pos[2]() + c->radius()*LJ_THRESH;
				lg_z = lg_z > 0 ? lg_z : 0;

				double le_x = c->pos[0]() - c->radius()*LJ_THRESH;
				if (le_x > LX) {
					le_x -= LX;
				}
				else if (le_x < 0) {
					le_x += LX;
				}
				double le_y = c->pos[1]() - c->radius()*LJ_THRESH;
				if (le_y > LY) {
					le_y -= LY;
				}
				else if (le_y < 0) {
					le_y += LY;
				}
				double le_z = c->pos[2]() - c->radius()*LJ_THRESH;
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
				if (lg_m_order == 0) {
					printf("eraaaa");
				}
				for (int i = DIV - 1; i >= 0; --i) {
					if (i < 0)printf("uyayuauh");
					if (((((unsigned)0b111) << (i * 3))&m_flag) != 0) {
						idx_in_shared = lg_m_order >> ((i + 1) * 3);
						break;
					}
					shared_level++;
				}
				
			
				//if (i > 98) {
					//printf("teuays");
				//}
				
				slist[((int)std::pow(8, shared_level) - 1) / 7 + idx_in_shared].push_back(c.get());
				
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