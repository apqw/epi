#include "cell_connection.h"
#include "CellManager.h"
#include "../../global.h"
#include "../../utils.h"
#include "../../misc/DynArr.h"
#include <iostream>
static void grid_init(CellManager& cman, Dyn3DArr<std::atomic<int>>& aindx, Dyn4DArr<Cell*>& area) { //DO NOT RESTRICT ptrs in area
    aindx.all_memset(0);

    cman.all_foreach_parallel_native([&](Cell*const c) {//cannot restrict
                                                        //auto&c = cman[i];
        int aix, aiy, aiz;
        aix = (int)((real(0.5)*pm->LX - p_diff_x(real(0.5)*pm->LX, c->x())) / pm->AREA_GRID);
        aiy = (int)((real(0.5)*pm->LY - p_diff_y(real(0.5)*pm->LY, c->y())) / pm->AREA_GRID);
        aiz = (int)((min0(c->z())) / pm->AREA_GRID);

        if ((aix >= (int)pm->ANX || aiy >= (int)pm->ANY || aiz >= (int)pm->ANZ || aix < 0 || aiy < 0 || aiz < 0)) {
            throw std::logic_error(
                "Bad cell position(on grid)\nRaw cell position: x=" + std::to_string(c->x()) + " y=" + std::to_string(c->y()) + " z=" + std::to_string(c->z()) + "\n"
                + "Grid cell position: x=" + std::to_string(aix) + " y=" + std::to_string(aiy) + " z=" + std::to_string(aiz) + "\n"
                + "Grid size: X=" + std::to_string(pm->ANX) + " Y=" + std::to_string(pm->ANY) + " Z=" + std::to_string(pm->ANZ) + "\n"
                + "Grid unit size:" + std::to_string(pm->AREA_GRID));
        }

        area.at(aix,aiy,aiz,aindx.at(aix,aiy,aiz)++) = c;
        c->connected_cell.force_set_count(c->state == MEMB ? pm->MEMB_ADJ_CONN_NUM : 0);
    });
}


static void connect_proc(CellManager& cman, const Dyn3DArr<std::atomic<int>>& aindx, const Dyn4DArr<Cell*>& area) {
#define LAT_ARR(axis)\
struct lat_arr_##axis {\
	std::vector<int> idx;\
	lat_arr_##axis():idx(pm->AN##axis * 3){\
        for (size_t i = 0; i < pm->AN##axis * 3; i++) {\
			idx[i] = i%pm-> AN##axis;\
		}\
	}\
	int operator[](int i) {\
		return idx[i];\
	}\
}
    LAT_ARR(X); LAT_ARR(Y); LAT_ARR(Z);
    static lat_arr_X xidx;
    static lat_arr_Y yidx;
    static lat_arr_Z zidx;

    constexpr int ii = 2;
        cman.non_memb_foreach_parallel_native([&](Cell*const c) { //cannot restrict

            const int anx = (int)(c->x() / pm->AREA_GRID);
            const int any = (int)(c->y() / pm->AREA_GRID);
            const int anz = (int)(c->z() / pm->AREA_GRID);

            if ((anx >= (int)pm->ANX || any >= (int)pm->ANY || anz >= (int)pm->ANZ || anx < 0 || any < 0 || anz < 0)) {
                throw std::logic_error(
                    "Bad cell position(on grid)\nRaw cell position: x=" + std::to_string(c->x()) + " y=" + std::to_string(c->y()) + " z=" + std::to_string(c->z()) + "\n"
                    + "Grid cell position: x=" + std::to_string(anx) + " y=" + std::to_string(any) + " z=" + std::to_string(anz)+"\n"
                    + "Grid size: X=" + std::to_string(pm->ANX) + " Y=" + std::to_string(pm->ANY) + " Z=" + std::to_string(pm->ANZ)+"\n"
                +"Grid unit size:"+std::to_string(pm->AREA_GRID));
            }
            const size_t my_index = c->get_index();
            const int xend = anx + ii + (int)pm->ANX;
            const int yend = any + ii + (int)pm->ANY;
            const int zend = anz + ii + (int)pm->ANZ;
            for (int j = anx - ii + (int)pm->ANX; j <= xend; j++) {
                const int aix = xidx[j];
                for (int k = any - ii + (int)pm->ANY; k <= yend; k++) {
                    const int aiy = yidx[k];
                    for (int l = anz - ii + (int)pm->ANZ; l <= zend; l++) {
                        const int aiz = zidx[l];
                        const int sz = aindx.at(aix, aiy, aiz);

                        for (int m = 0; m < sz; ++m) {
                            Cell*const o = area.at(aix, aiy, aiz, m);
                            if (my_index <= o->get_index())continue;
                            const real diffx = p_diff_x(c->x(), o->x());
                            const real diffy = p_diff_y(c->y(), o->y());
                            const real diffz = c->z() - o->z();
                            const real rad_sum = c->radius + o->radius;
                            if (DIST_SQ(diffx, diffy, diffz) <= pm->LJ_THRESH*pm->LJ_THRESH*rad_sum*rad_sum) {
                                c->connected_cell.push_back(o);
                                o->connected_cell.push_back(c);
                                //assert(c->connected_cell.size() < N2);
                            }
                        }
                    }
                }
            }
        });
}

static const Cell* find_dermis(const Cell*const RESTRICT c) {
    real d1Sq = pm->LX*pm->LX;
    real distanceSq = 0;
    const Cell* dermis = nullptr;
    c->connected_cell.foreach([&](const Cell*const cptr) { //no restrict
        if (cptr->state == MEMB) {
            distanceSq = p_cell_dist_sq(c, cptr);
            if (distanceSq < d1Sq) {//2��ɂȂ��ĂȂ�����
                d1Sq = distanceSq;
                dermis = cptr;
            }
        }
    });
    return dermis;
}

static void set_dermis(CellManager& cman) {
    cman.other_foreach_parallel_native([](Cell*const RESTRICT c) {

        if (c->state == FIX || c->state == MUSUME) {
            c->set_dermis(find_dermis(c));
        }
    });
}

void connect_cell(CellManager& cman) {
    static Dyn3DArr<std::atomic<int>> aindx(pm->ANX, pm->ANY, pm->ANZ);
    static Dyn4DArr<Cell*> area(pm->ANX, pm->ANY, pm->ANZ,pm->N3);
    grid_init(cman, aindx, area);
    connect_proc(cman, aindx, area);
    set_dermis(cman);
}
