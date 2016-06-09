#include "cell_connection.h"
#include "define.h"
#include "utils.h"
#include "cell.h"
#include "cell_conn_value.h"
#include <cinttypes>
using cint = int_fast16_t;
#define CINT_F SCNiFAST16

void grid_init(CellManager& cman, std::atomic<cint>(&aindx)[cont::ANX][cont::ANY][cont::ANZ] , Cell* RESTRICT(&area)[cont::ANX][cont::ANY][cont::ANZ][cont::N3]) {
	using namespace cont;
	std::memset(aindx, 0, sizeof(std::atomic<cint>)*ANX*ANY*ANZ);
	cman.all_foreach_parallel_native([&](Cell*const RESTRICT& c) {
		//auto&c = cman[i];
		cint aix, aiy, aiz;
		aix = (0.5*LX - p_diff_x(0.5*LX, c->x())) / AREA_GRID;
		aiy = (0.5*LY - p_diff_y(0.5*LY, c->y())) / AREA_GRID;
		aiz = (min0(c->z())) / AREA_GRID;

        if ((aix >= (cint)ANX || aiy >= (cint)ANY || aiz >= (cint)ANZ || aix < 0 || aiy < 0 || aiz < 0)) {
			printf("err\n");
			printf("cx:%lf cy:%lf cz:%lf\n", c->x(), c->y(), c->z());
            printf("aix:%" CINT_F " aiy:%" CINT_F " aiz:%" CINT_F "\n", aix, aiy, aiz);
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

void connect_proc(CellManager& cman, const std::atomic<cint>(&aindx)[cont::ANX][cont::ANY][cont::ANZ], Cell*const RESTRICT(&area)[cont::ANX][cont::ANY][cont::ANZ][cont::N3]) {
	
#define LAT_ARR(axis)\
struct lat_arr_##axis {\
	cint idx[cont::AN##axis * 3];\
	lat_arr_##axis() {\
        for (size_t i = 0; i < cont::AN##axis * 3; i++) {\
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
    constexpr cint ii = 2;
	cman.non_memb_foreach_parallel_native([&](Cell*const&c) { //cannot restrict
		
		const cint anx = (cint)(c->x() / AREA_GRID);
		const cint any = (cint)(c->y() / AREA_GRID);
		const cint anz = (cint)(c->z() / AREA_GRID);

        assert(!(anx >= (cint)ANX || any >= (cint)ANY || anz >= (cint)ANZ || anx < 0 || any < 0 || anz < 0));

        for (cint j = anx - ii+(cint)ANX; j <= anx + ii+(cint)ANX; j++) {
			const cint aix = xidx[j];
            for (cint k = any - ii + (cint)ANY; k <= any + ii + (cint)ANY; k++) {
				const cint aiy = yidx[k];
                for (cint l = anz - ii + (cint)ANZ; l <= anz + ii + (cint)ANZ; l++) {
					const cint aiz = zidx[l];
					const cint sz = aindx[aix][aiy][aiz];

					for (cint m = 0; m < sz; ++m) {
						Cell*const o = area[aix][aiy][aiz][m];
						if (c->get_index() <= o->get_index())continue;
						const double diffx = p_diff_x(c->x(), o->x());
						const double diffy = p_diff_y(c->y(), o->y());
						const double diffz = c->z() - o->z();
						const double rad_sum = c->radius + o->radius;
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

	cman.all_foreach_parallel_native([](Cell*const RESTRICT c) {

        for(auto&& it = c->gj.begin(),itend=c->gj.end();it!=itend;++it){
            it->second.uncheck();
        }

        c->connected_cell.foreach([&c](Cell*const RESTRICT cptr) {
            c->gj[cptr].check();
        });

        for(auto&& it = c->gj.begin(),itend=c->gj.end();it!=itend;++it){
            if(!it->second.is_checked()){
                it->second.reset();
            }
        }

        /*
		for (auto&& it = c->gj.begin(); it != c->gj.end();) {
			if (!c->connected_cell.exist(it->first)) {
				it = c->gj.erase(it);
			}
			else {
				++it;
			}
		}

		c->connected_cell.foreach([c](Cell* cptr) {
			if (c->gj.find(cptr)==c->gj.end()) {
				c->gj.emplace(cptr, cont::gj_init);
			}
		});
    */
	});
}

const Cell* find_dermis(const Cell*const RESTRICT c) {
	double d1Sq = cont::LX*cont::LX;
	double distanceSq = 0;
	const Cell* dermis = nullptr;
	c->connected_cell.foreach([&](const Cell*const RESTRICT& cptr) {
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

inline void set_dermis(CellManager& cman) {
	cman.other_foreach_parallel_native([](Cell*const RESTRICT&c) {
		
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
