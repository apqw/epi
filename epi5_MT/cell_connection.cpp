#include "cell_connection.h"
#include "define.h"
#include "utils.h"
#include "cell.h"

using cint = int_fast16_t;

void grid_init(CellManager& cman, std::atomic<cint>(&aindx)[cont::ANX][cont::ANY][cont::ANZ] , Cell*(&area)[cont::ANX][cont::ANY][cont::ANZ][cont::N3]) {
	using namespace cont;
	std::memset(aindx, 0, sizeof(std::atomic<cint>)*ANX*ANY*ANZ);
	cman.all_foreach_parallel([&](size_t i) {
		auto&c = cman[i];
		cint aix, aiy, aiz;
		aix = (0.5*LX - p_diff_x(0.5*LX, c->x())) / AREA_GRID;
		aiy = (0.5*LY - p_diff_y(0.5*LY, c->y())) / AREA_GRID;
		aiz = (min0(c->z())) / AREA_GRID;

		if ((aix >= ANX || aiy >= ANY || aiz >= ANZ || aix < 0 || aiy < 0 || aiz < 0)) {
			printf("err\n");
			printf("cx:%lf cy:%lf cz:%lf\n", c->x(), c->y(), c->z());
			printf("aix:%d aiy:%d aiz:%d\n", aix, aiy, aiz);
			assert(false);
		}

		area[aix][aiy][aiz][aindx[aix][aiy][aiz]++] = c;
		/*
		Ç±ÇÃassertÇÕè¡ÇµÇƒÇÊÇ¢
		*/
		assert(aindx[aix][aiy][aiz] < N3);
		c->connected_cell.force_set_count(c->state == MEMB ? 4 : 0);

	});
}

void connect_proc(CellManager& cman, std::atomic<cint>(&aindx)[cont::ANX][cont::ANY][cont::ANZ], Cell*(&area)[cont::ANX][cont::ANY][cont::ANZ][cont::N3]) {
	
#define LAT_ARR(axis)\
struct lat_arr_##axis {\
	cint idx[cont::AN##axis * 3];\
	lat_arr_##axis() {\
		for (int i = 0; i < cont::AN##axis * 3; i++) {\
			idx[i] = i%cont::AN##axis;\
		}\
	}\
	cint operator[](int i) {\
		return idx[i];\
	}\
}
	LAT_ARR(X); LAT_ARR(Y); LAT_ARR(Z);

	static lat_arr_X xidx;
	static lat_arr_Y yidx;
	static lat_arr_Z zidx;
	

	using namespace cont;
	constexpr int ii = 2;
	cman.non_memb_foreach_parallel([&](size_t i) {
		auto&c = cman[i];
		cint anx = (cint)(c->x() / AREA_GRID);
		cint any = (cint)(c->y() / AREA_GRID);
		cint anz = (cint)(c->z() / AREA_GRID);

		assert(!(anx >= ANX || any >= ANY || anz >= ANZ || anx < 0 || any < 0 || anz < 0));

		for (int j = anx - ii+ANX; j <= anx + ii+ANX; j++) {
			cint aix = xidx[j];
			for (int k = any - ii + ANY; k <= any + ii + ANY; k++) {
				cint aiy = yidx[k];
				for (int l = anz - ii + ANZ; l <= anz + ii + ANZ; l++) {
					cint aiz = zidx[l];
					cint sz = aindx[aix][aiy][aiz];

					for (cint m = 0; m < sz; ++m) {
						Cell* o = area[aix][aiy][aiz][m];
						if (c->get_index() <= o->get_index())continue;
						double diffx = p_diff_x(c->x(), o->x());
						double diffy = p_diff_y(c->y(), o->y());
						double diffz = c->z() - o->z();
						double rad_sum = c->radius + o->radius;
						if (DIST_SQ(diffx, diffy, diffz) <= LJ_THRESH*LJ_THRESH*rad_sum*rad_sum) {
							c->connected_cell.push_back(o);
							o->connected_cell.push_back(c);
							assert(c->connected_cell.size() < N2);
						}
					}
				}
			}
		}
	});
}

void gj_refresh(CellManager& cman) {
	static std::unordered_map<Cell*,bool> connflg[cont::MAX_CELL_NUM];

	cman.all_foreach_parallel([&](size_t i) {
		auto&c = cman[i];
		
		auto&cf = connflg[i];
		for (auto&& it = cf.begin(); it != cf.end(); ++it) {
			it->second = false;
		}
		c->connected_cell.foreach([&](Cell* cptr) {
			cf[cptr] = true;
			if (!c->gj.exist(cptr)) {
				c->gj._emplace(cptr, cont::gj_init);
			}
		});

		for (auto&& it = cf.begin(); it != cf.end(); ++it) {
			if (!it->second) {
				c->gj._set(it->first, cont::gj_init);
			}
		}
	
	});
}

Cell* find_dermis(Cell* c) {
	double d1Sq = cont::LX*cont::LX;
	double distanceSq = 0;
	Cell* dermis = nullptr;
	c->connected_cell.foreach([&](Cell* cptr) {
		if (cptr->state == MEMB) {
			distanceSq = p_cell_dist_sq(c, cptr);
			if (distanceSq < d1Sq) {//2èÊÇ…Ç»Ç¡ÇƒÇ»Ç©Ç¡ÇΩ
				d1Sq = distanceSq;
				dermis = cptr;
			}
		}
	});
	return dermis;
}

void set_dermis(CellManager& cman) {
	cman.other_foreach_parallel([&](size_t i) {
		auto&c = cman[i];
		if (c->state == FIX || c->state == MUSUME) {
			c->set_dermis(find_dermis(c));
		}
	});
}

void connect_cell(CellManager & cman)
{
	using namespace cont;
	static std::atomic<cint> aindx[ANX][ANY][ANZ] = {};
	static Cell* area[ANX][ANY][ANZ][N3] = { nullptr };
	

	grid_init(cman,aindx,area);
	connect_proc(cman, aindx, area);
	gj_refresh(cman);
	set_dermis(cman);
	
}
