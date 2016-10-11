#pragma once
#include "../../global.h"
#include "../../define.h"
#include "atomics.h"
#include "DualValue.h"
#include <string>
#include "../../util/vec/Vec.h"
class Cell
{
private:
    const Cell* _dermis=nullptr;
    size_t index;
    struct ctor_cookie{};

public:
    friend class CellManager;
    struct Memb_data{
        	/*
        	 * counterclockwise (starts from +x axis direction)
        	 */
        	std::vector<Cell*> memb;

        	real nv[3]; real ipn;
        	real mv[3]; real ipm;
        	real dn, dm;

        	real nv_a[3]; real ipn_a;
        	real mv_a[3]; real ipm_a;
        	real dn_a, dm_a;
        } md;
    dual_real x, y, z;
    Lockfree_push_stack_dyn<Cell*> connected_cell;
    CELL_STATE state;
    const int fix_origin;
    Cell*pair = nullptr;
    real ca2p;
    real ca2p_avg;
    real IP3;
    real ex_inert;
    real agek;
    real ageb;
    real ex_fat;
    real in_fat;
    real spr_nat_len;
    real radius;
    unsigned int rest_div_times;
    bool is_malignant;
    Vec<3,int> lat;
    real diff_u;
    real div_age_thresh;
    bool is_touch;
    real spring_force_to_memb;
    const Cell* dermis()const;
    void set_dermis(const Cell *const);
    static real get_div_age_thresh(CELL_STATE state);
    static int correct_div_times(CELL_STATE state, int given_times);
    Cell(ctor_cookie,CELL_STATE _state, int _fix_origin,
        real _ca2p = pm->ca2p_init,
        real _IP3 = pm->IP3_init,
        real _ex_inert = 0,
        real _agek = 0, real _ageb = 0, real _ex_fat = 0, real _in_fat = 0, real _spr_nat_len = 0,
        real _x = 0, real _y = 0, real _z = 0, int _rest_div_time = 0,
        real _radius = pm->R_max, real _ca2p_avg = pm->ca2p_init,
        bool _is_malignant = false);
    void set_index(size_t idx);
    size_t get_index()const;
    void migrate(size_t destidx);
    std::string cell_info_str();
    Vec<3, real> pos_as_vector()const;
};

/** ダブルカウントの回避条件 */
inline bool no_double_count(const Cell*const c1, const Cell*const c2) {
    return c1->get_index() > c2->get_index();
}

inline real p_cell_dist_sq(const Cell*const RESTRICT c1, const Cell*const RESTRICT c2)
{
    return p_dist_sq(c1->x(), c1->y(), c1->z(), c2->x(), c2->y(), c2->z());
}
