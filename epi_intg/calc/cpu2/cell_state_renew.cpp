
#define _USE_MATH_DEFINES
#include <cmath>
#include "cell_state_renew.h"
#include "CellManager.h"
#include "../../global.h"
#include "../../utils.h"
#include "../../util/rand/Random_gen.h"
#include "../../util/vec/Vec.h"

////////////////////////////////////////////////////////////////////////////////
//
//    agebの計算の定義
//
////////////////////////////////////////////////////////////////////////////////

/**
 *  @param [in] c 計算対象の細胞細胞
 *  @return 条件によって補正されたeps_kb
 */
inline real weighted_eps_kb(const Cell*const c) {
    return c->is_malignant ? pm->accel_div *pm->eps_kb : pm->eps_kb;
}

inline real ageb_const(const Cell*const RESTRICT c) {

    return weighted_eps_kb(c)*(pm->S2 + pm->alpha_b*min0(c->ca2p_avg - pm->ca2p_init));
}



////////////////////////////////////////////////////////////////////////////////
//
//    agekの計算の定義
//
////////////////////////////////////////////////////////////////////////////////

/**
 *  @param [in] c 計算対象の細胞細胞
 *  @return 条件によって補正されたeps_ks
 */
inline real weighted_eps_ks(const Cell*const c) {
    return c->is_malignant ? pm->accel_diff *pm->eps_ks : pm->eps_ks;
}

inline real agek_ALIVE_const(const Cell*const RESTRICT c) {

    return weighted_eps_ks(c)*(pm->S0 + pm->alpha_k*min0(c->ca2p_avg - pm->ca2p_init));
}

inline real agek_DEAD_AIR_const() {

    return pm->eps_kk*pm->S1;
}

////////////////////////////////////////////////////////////////////////////////
//
//    脂質の計算の定義
//
////////////////////////////////////////////////////////////////////////////////

/**
 *  @param [in] c 計算対象の細胞細胞
 *  @return 放出する脂質の時間差分
 */
inline real k_lipid_release(const Cell*const RESTRICT c) {
    return 0.25*pm->lipid_rel*(1 + tanh((c->ca2p_avg - pm->ubar) / pm->delta_sig_r1))*(1 + tanh((c->agek - pm->THRESH_SP) / pm->delta_lipid));
}

/**
 *  @param [in] c 計算対象の細胞細胞
 *  @return 生成する脂質の時間差分
 *  @attention この値が直接生成量になるわけではない。
 */
inline real k_lipid(const Cell*const RESTRICT c) {
    return 0.25*pm->lipid*(1 + tanh((pm->ubar - c->ca2p_avg) / pm->delta_sig_r1))*(1 + tanh((c->agek - pm->THRESH_SP) / pm->delta_lipid));
}




////////////////////////////////////////////////////////////////////////////////
//
//    細胞分裂に関する定義
//
////////////////////////////////////////////////////////////////////////////////

/**
 *  分裂するかどうかを確率的に決定する部分
 *  @return 分裂を行う場合true、そうでないならfalse
 *  @attention 確率的な分裂を行わない場合、常にtrueを返す。
 */
inline bool stochastic_div_test(const Cell*const c) {
    if (!pm->STOCHASTIC) {
        return true;
    }
    else {
        return genrand_real()*(c->div_age_thresh*pm->stoch_div_time_ratio) <= pm->DT_Cell*weighted_eps_kb(c)*pm->S2;
    }
}

/**
    分裂の準備ができているかどうか
    @attention 確率的な分裂を行う場合、乱数により返り値は変わる。
*/
inline bool is_divide_ready(const Cell*const RESTRICT c) {

    if (c->pair == nullptr && (c->ageb >= c->div_age_thresh*(1.0 - pm->stoch_div_time_ratio))) {
        return stochastic_div_test(c);
    }
    else {
        return false;
    }
}

inline Vec<3,real> calc_cell_uvec(const Cell*const RESTRICT c1, const Cell*const RESTRICT c2) {

	return vpsub(c1->pos_as_vector(),c2->pos_as_vector()).normalize();

}
Vec<3,real> div_direction(const Cell*const RESTRICT me, const Cell*const RESTRICT dermis) {

    Vec<3,real> nv=calc_cell_uvec(me, dermis);
    Vec<3,real> ov;
    real sum;
    do {
        real rand_theta = M_PI*genrand_real();
        real cr1 = cos(rand_theta);
        real sr1 = sin(rand_theta);
        real rand_phi = 2 * M_PI*genrand_real();
        real cr2 = cos(rand_phi);
        real sr2 = sin(rand_phi);
        Vec<3,real> rv ={sr1*cr2,sr1*sr2,cr1};
ov=Vec<3,real>::cross(nv,rv);
    } while ((sum = ov.norm_sq()) < 1.0e-14);
    return ov/sqrt(sum);
}

void cell_divide(CellManager& cman, Cell*const RESTRICT div) {

    if (div->dermis() == nullptr) {
    	throw std::runtime_error("No dermis found in divide_try."_s+div->cell_info_str());
    }

    div->pair = cman.create(
        MUSUME,
        div->fix_origin,
        div->x(), div->y(), div->z(),
        div->radius,
        div->ca2p,
        div->ca2p,//this is avg value,do not use orignal avg
        div->IP3,
        div->ex_inert,
        0, 0,//set ages 0
        0, 0,//set fats 0
        pm->delta_L,//init
        div->rest_div_times,
        div->is_malignant
        );
    div->pair->pair = div;
    div->spr_nat_len = pm->delta_L;
    div->ageb = 0;
    if (div->state == MUSUME) {
        div->rest_div_times--;
        div->pair->rest_div_times--;
    }

    Vec<3,real>dv=div_direction(div, div->dermis())*(0.5*pm->delta_L);
    //!set value
    div->x._set(div->x() + dv[0]);
    div->y._set(div->y() + dv[1]);
    div->z._set(div->z() + dv[2]);

    div->pair->x._set(div->pair->x() - dv[0]);
    div->pair->y._set(div->pair->y() - dv[1]);
    div->pair->z._set(div->pair->z() - dv[2]);
    std::cout<<"New cell:"+std::to_string(cman.size())<<std::endl;
}


////////////////////////////////////////////////////////////////////////////////
//
//    細胞の状態更新の定義
//
////////////////////////////////////////////////////////////////////////////////

/**
 *  娘細胞の状態更新。
 *  分化、分裂、細胞周期の更新を行う。
 */
void _MUSUME_state_renew(CellManager& cman, Cell*const RESTRICT musume) {
    if (musume->dermis() == nullptr && musume->pair == nullptr) {
        if (pm->SYSTEM == WHOLE) {
            std::cout<<"ALIVE detected"<<std::endl;
            musume->state = ALIVE;
            //musume->cst.k_aging_start_timestep=cman.current_timestep;
        }
        else if (pm->SYSTEM == BASAL) {
            //musume->state = DISA;
            cman.add_remove_queue(musume->get_index());
        }
        return;
    }

    if (musume->rest_div_times > 0 && is_divide_ready(musume)) {
        cell_divide(cman, musume);
    }
    else {
        musume->ageb += pm->DT_Cell*ageb_const(musume);
    }
}

/**
 *  幹細胞の状態更新。
 *  分裂、細胞周期の更新を行う。
 */
void _FIX_state_renew(CellManager& cman, Cell*const RESTRICT fix) {
    if (fix->dermis() == nullptr) {
        throw std::runtime_error("No dermis found in FIX_state_renew."_s+fix->cell_info_str());
    }

    if (is_divide_ready(fix)) {
        cell_divide(cman, fix);
    }
    else {
        fix->ageb += pm->DT_Cell*ageb_const(fix);
    }
}

/**
 *  角層、空気の状態更新。
 *  加齢、剥離を行う。
 */
void _DEAD_AIR_state_renew(CellManager& cman, Cell*const RESTRICT da) {

    if (da->agek >= pm->ADHE_CONST&&da->connected_cell.size() <= pm->DISA_conn_num_thresh) {
        //da->state = DISA;
        cman.add_remove_queue(da->get_index());
    }
    else {
        da->agek += pm->DT_Cell*agek_DEAD_AIR_const();
    }
}

/**
 *  角化を行う。
 *  ALIVE->DEAD
 */
void cornificate(CellManager & cman, Cell * const RESTRICT al)
{
	al->state = DEAD;
   // al->cst.k_cornified_timestep=cman.current_timestep;
	std::cout<<"sw updated"<<++cman.sw<<std::endl;
}
/**
 *  有棘細胞の状態更新。
 *  加齢、角化を行う。
 */
inline void _ALIVE_state_renew(CellManager& cman, Cell*const RESTRICT al) {
    if (al->agek >= pm->THRESH_DEAD) {
        cornificate(cman, al);
    }
    else {
        const real tmp = k_lipid_release(al)*al->in_fat;
        al->in_fat += pm->DT_Cell*(k_lipid(al)*(1.0 - al->in_fat) - tmp);
        al->ex_fat += pm->DT_Cell*tmp;
        al->agek += pm->DT_Cell*agek_ALIVE_const(al); //update last
    }
}

/**
 *  分裂中の細胞の状態更新。
 *  ばねの自然長の更新、分裂の完了処理を行う。
 */
void pair_disperse(Cell*const c) {//cannot restrict due to c->pair->pair (== c)
if(c->pair==nullptr){
	throw std::logic_error("Pair cell does not exist."_s + c->cell_info_str());
}

if(c->pair->pair!=c){
	throw std::logic_error("Inconsistent pairing."_s + c->cell_info_str());
}
    real rad_sum = c->radius + c->pair->radius;
    real unpair_th = pm->unpair_dist_coef*rad_sum;
    real distSq = 0;
    if (c->spr_nat_len < 2.0*c->radius) {
        c->spr_nat_len += pm->DT_Cell*pm->eps_L;
        c->pair->spr_nat_len = c->spr_nat_len;
    }
    else if ((distSq = p_cell_dist_sq(c, c->pair))>unpair_th*unpair_th) {
        c->spr_nat_len = 0;
        c->pair->spr_nat_len = 0;
        c->pair->pair = nullptr;
        c->pair = nullptr;
        std::cout<<"Unpaired."<<std::endl;
    }
}


////////////////////////////////////////////////////////////////////////////////
//
//    外部公開関数
//
////////////////////////////////////////////////////////////////////////////////


/**
 *  全ての細胞の状態を更新する。
 *  @attention 並列化されている。
 */
void cell_state_renew(CellManager & cman)
{
    cman.other_foreach_parallel_native([&](Cell*const RESTRICT c) {
        //auto&c = cman[i];
        switch (c->state) {
        case ALIVE:
            _ALIVE_state_renew(cman, c);
            break;
        case DEAD:case AIR:
            _DEAD_AIR_state_renew(cman, c);
            break;
        case FIX:
            _FIX_state_renew(cman, c);
            break;
        case MUSUME:
            _MUSUME_state_renew(cman, c);
            break;
        default:
            break;
        }
    });
    cman.remove_exec();
    cman.other_foreach_parallel_native([](Cell*const RESTRICT c) {

        if (c->pair != nullptr) if (no_double_count(c, c->pair)) {
                pair_disperse(c);
        }

    });
}

/**
 *  強制的にカルシウム刺激を与える場合の状態の初期化を行う。
 *  @param cman CellManager
 *  @param [in] zzmax 細胞が存在する最大高さ(z座標)
 */
void initialize_sc(CellManager & cman, real zzmax)
{
    cman.other_foreach_parallel_native([zzmax](Cell*const RESTRICT c) {
        if (c->state == ALIVE&&zzmax - c->z() < 8.0 * pm->R_max) {
            c->agek = c->agek > pm->THRESH_DEAD ? c->agek : pm->THRESH_DEAD;
            c->state = AIR;
        }
    });
}
