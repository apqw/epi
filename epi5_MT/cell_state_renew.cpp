#include "cell_state_renew.h"
#include "cell.h"
#include "cellmanager.h"
#include "define.h"
#include "Random_gen.h"
#include "utils.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <cassert>
#include <iostream>
/**
 *  @file 細胞の状態更新に関する操作を定義。
 */


static constexpr double stoch_div_time_ratio = 0.25;

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//    agebの計算の定義                                                        //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

static constexpr double S2 = 0.1; //TODO:よりわかりやすい命名

/**
 *  @param [in] c 計算対象の細胞細胞
 *  @return 条件によって補正されたeps_kb
 */
inline double weighted_eps_kb(const Cell*const c) {
    static constexpr double accel_div = 1.0;
    static constexpr double eps_kb = 0.12;
    return c->is_malignant ? accel_div *eps_kb : eps_kb;
}

inline double ageb_const(const Cell*const RESTRICT c) {
    using namespace cont;
    static constexpr double alpha_b = 5.0;

    return weighted_eps_kb(c)*(S2 + alpha_b*min0(c->ca2p_avg - ca2p_init));
}


////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//    agekの計算の定義                                                        //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

static constexpr double S0 = 0.1*0.2;//TODO:よりわかりやすい命名

/**
 *  @param [in] c 計算対象の細胞細胞
 *  @return 条件によって補正されたeps_ks
 */
inline double weighted_eps_ks(const Cell*const c) {
    static constexpr double eps_ks = 0.10*0.5;
    static constexpr double accel_diff = 1.0;
    return c->is_malignant ? accel_diff *eps_ks : eps_ks;
}

inline double agek_ALIVE_const(const Cell*const RESTRICT c) {
    using namespace cont;
    static constexpr double alpha_k = 2.0;

    return weighted_eps_ks(c)*(S0 + alpha_k*min0(c->ca2p_avg - ca2p_init));
}

inline double agek_DEAD_AIR_const() {
    static constexpr double eps_kk = 0.10*0.5;
    static constexpr double S1 = 0.1;

    return eps_kk*S1;
}

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//    脂質の計算の定義                                                        //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

/** 脂質の生成・放出を切り替えるカルシウム濃度 */
static constexpr double ubar = 0.45;

/**
 *  agekによるスイッチングの緩さ
 *  @note tanhの分母
 */
static constexpr double delta_lipid = 0.1;

/**
 *  カルシウム濃度によるスイッチングの緩さ
 *  @note tanhの分母
 */
static constexpr double delta_sig_r1 = 0.1;

/**
 *  @param [in] c 計算対象の細胞細胞
 *  @return 放出する脂質の時間差分
 */
inline double k_lipid_release(const Cell*const RESTRICT c) {
    using namespace cont;

    static constexpr double lipid_rel = 0.05*4.0;
    return 0.25*lipid_rel*(1 + tanh((c->ca2p_avg - ubar) / delta_sig_r1))*(1 + tanh((c->agek - THRESH_SP) / delta_lipid));
}

/**
 *  @param [in] c 計算対象の細胞細胞
 *  @return 生成する脂質の時間差分
 *  @attention この値が直接生成量になるわけではない。
 */
inline double k_lipid(const Cell*const RESTRICT c) {
    using namespace cont;

    static constexpr double lipid = 0.6;
    return 0.25*lipid*(1 + tanh((ubar - c->ca2p_avg) / delta_sig_r1))*(1 + tanh((c->agek - THRESH_SP) / delta_lipid));
}


////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//    細胞分裂に関する定義                                                    //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

/**
 *  分裂するかどうかを確率的に決定する部分
 *  @return 分裂を行う場合true、そうでないならfalse
 *  @attention 確率的な分裂を行わない場合、常にtrueを返す。
 */
inline bool stochastic_div_test(const Cell*const c) {
    using namespace cont;
    if (!STOCHASTIC) {
        return true;
    }
    else {
        return genrand_real()*(c->div_age_thresh*stoch_div_time_ratio) <= stoch_corr_coef*DT_Cell*weighted_eps_kb(c)*S2;
    }
}

/**
    分裂の準備ができているかどうか
    @attention 確率的な分裂を行う場合、乱数により返り値は変わる。
*/
inline bool is_divide_ready(const Cell*const RESTRICT c) {

    if (c->pair == nullptr && (c->ageb >= c->div_age_thresh*(1.0 - stoch_div_time_ratio))) {
        return stochastic_div_test(c);
    }
    else {
        return false;
    }
}

/**
 *  c2からc1への単位ベクトルを計算する。
 *  @param [in] c1 細胞c1
 * 	@param [in] c2 細胞c2
 *  @param [out] outx 単位ベクトルのx成分
 *  @param [out] outy 単位ベクトルのy成分
 *  @param [out] outz 単位ベクトルのz成分
 *
 *  @attention c1とc2が同じ細胞を指している場合、動作は未定義。
 *  また、outx,outy,outzのいずれかが同じ変数を指している場合も動作は未定義。
 */
inline void calc_cell_uvec(const Cell*const RESTRICT c1, const Cell*const RESTRICT c2, double*const RESTRICT outx, double*const RESTRICT outy, double*const RESTRICT outz) {
    double nvx = p_diff_x(c1->x(), c2->x());
    double nvy = p_diff_y(c1->y(), c2->y());
    double nvz = c1->z() - c2->z();

    double norm = sqrt(DIST_SQ(nvx, nvy, nvz));
    *outx = nvx / norm;
    *outy = nvy / norm;
    *outz = nvz / norm;

}

/**
 *  最も近いdermisに対しての分裂方向(単位ベクトル)をランダムに計算する。
 *
 *  @param [in] me 分裂を行う細胞
 *  @param [in] dermis 最近傍のdermis
 *  @param [out] outx 単位ベクトルのx成分
 *  @param [out] outy 単位ベクトルのy成分
 *  @param [out] outz 単位ベクトルのz成分
 *
 *  @attention meとdermisが同じ細胞を指している場合、動作は未定義。
 *  また、outx,outy,outzのいずれかが同じ変数を指している場合も動作は未定義。
 */
void div_direction(const Cell*const RESTRICT me, const Cell*const RESTRICT dermis, double*const RESTRICT outx, double*const RESTRICT outy, double*const RESTRICT outz) {

    double nx, ny, nz;
    calc_cell_uvec(me, dermis, &nx, &ny, &nz);
    double ox, oy, oz;
    double sum;
    do {
        double rand_theta = M_PI*genrand_real();
        double cr1 = cos(rand_theta);
        double sr1 = sin(rand_theta);
        double rand_phi = 2 * M_PI*genrand_real();
        double cr2 = cos(rand_phi);
        double sr2 = sin(rand_phi);
        double rvx = sr1*cr2;
        double rvy = sr1*sr2;
        double rvz = cr1;

        ox = ny*rvz - nz*rvy;
        oy = nz*rvx - nx*rvz;
        oz = nx*rvy - ny*rvx;
    } while ((sum = DIST_SQ(ox, oy, oz)) < 1.0e-14);
    sum = sqrt(sum);
    *outx = ox / sum;
    *outy = oy / sum;
    *outz = oz / sum;
}

/**
 *  細胞分裂を行う。
 *
 *  @param cman cellmanager
 *  @param div 分裂を行う細胞
 *
 *  @attention 最近傍のdermisが計算されている必要がある。
 */
void cell_divide(CellManager& cman, Cell*const RESTRICT div) {


    using namespace cont;
    static constexpr double delta_L = 0.01*cont::R_max;
    div->pair = cman.create(
        MUSUME,
        div->fix_origin,
        div->x(), div->y(), div->z(),
        div->radius,
        div->ca2p(),
        div->ca2p(),//this is avg value,do not use orignal avg
        div->IP3(),
        div->ex_inert,
        0, 0,//set ages 0
        0, 0,//set fats 0
        delta_L,//init
        div->rest_div_times,
        div->is_malignant
        );
    div->pair->pair = div;
    div->spr_nat_len = delta_L;
    div->ageb = 0;
    if (div->state == MUSUME) {
        div->rest_div_times--;
        div->pair->rest_div_times--;
    }
    if (div->dermis() == nullptr) {
        std::cout << "No dermis found in divide_try." << std::endl;
        exit(1);
    }
    //assert(div->dermis() != nullptr);
    double divx, divy, divz;
    div_direction(div, div->dermis(), &divx, &divy, &divz);
    //!set value
    div->x._set(div->x() + divx*0.5*delta_L);
    div->y._set(div->y() + divy*0.5*delta_L);
    div->z._set(div->z() + divz*0.5*delta_L);

    div->pair->x._set(div->pair->x() - divx*0.5*delta_L);
    div->pair->y._set(div->pair->y() - divy*0.5*delta_L);
    div->pair->z._set(div->pair->z() - divz*0.5*delta_L);
    printf("new cell detected. cell_num:%zd\n",cman.size());
}

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//    細胞の状態更新の定義                                                    //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

/**
 *  娘細胞の状態更新。
 *  分化、分裂、細胞周期の更新を行う。
 */
void _MUSUME_state_renew(CellManager& cman, Cell*const RESTRICT musume) {
    if (musume->dermis() == nullptr && musume->pair == nullptr) {
        if (SYSTEM == WHOLE) {
            printf("ALIVE detected\n");
            musume->state = ALIVE;
        }
        else if (SYSTEM == BASAL) {
            musume->state = DISA;
            cman.add_remove_queue(musume->get_index());
        }
        return;
    }

    if (musume->rest_div_times > 0 && is_divide_ready(musume)) {
        cell_divide(cman, musume);
    }
    else {
        musume->ageb += cont::DT_Cell*ageb_const(musume);
    }
}

/**
 *  幹細胞の状態更新。
 *  分裂、細胞周期の更新を行う。
 */
void _FIX_state_renew(CellManager& cman, Cell*const RESTRICT fix) {
    if (fix->dermis() == nullptr) {
        printf("err\n");
        printf("x:%lf,y:%lf,z:%lf\n", fix->x(), fix->y(), fix->z());
        printf("connected_num:%zd\n", fix->connected_cell.size());
        assert(fix->dermis() != nullptr);
        exit(1);
        return;
    }

    if (is_divide_ready(fix)) {
        cell_divide(cman, fix);
    }
    else {
        fix->ageb += cont::DT_Cell*ageb_const(fix);
    }
}

/**
 *  角層、空気の状態更新。
 *  加齢、剥離を行う。
 */
void _DEAD_AIR_state_renew(CellManager& cman, Cell*const RESTRICT da) {
    static constexpr double ADHE_CONST = 31.3;
    static constexpr unsigned int DISA_conn_num_thresh = 11; //Nc

    if (da->agek >= ADHE_CONST&&da->connected_cell.size() <= DISA_conn_num_thresh) {
        da->state = DISA;
        cman.add_remove_queue(da->get_index());
    }
    else {
        da->agek += cont::DT_Cell*agek_DEAD_AIR_const();
    }
}

/**
 *  有棘細胞の状態更新。
 *  加齢、角化を行う。
 */
inline void _ALIVE_state_renew(CellManager& cman, Cell*const RESTRICT al) {
    using namespace cont;
    if (al->agek >= THRESH_DEAD) {
        cornificate(cman, al);
    }
    else {
        const double tmp = k_lipid_release(al)*al->in_fat;
        al->in_fat += DT_Cell*(k_lipid(al)*(1.0 - al->in_fat) - tmp);
        al->ex_fat += DT_Cell*tmp;
        al->agek += DT_Cell*agek_ALIVE_const(al); //update last
    }
}

/**
 *  分裂中の細胞の状態更新。
 *  ばねの自然長の更新、分裂の完了処理を行う。
 */
void pair_disperse(Cell*const c) {//cannot restrict due to c->pair->pair (== c)

    static constexpr double eps_L = 0.14;//ok
    static constexpr double unpair_dist_coef = 0.9;
    //fprintf(stdout, "disperse start.\n");
    //fflush(stdout);
    using namespace cont;
    assert(c->pair != nullptr);
    assert(c->pair->pair == c);
    double rad_sum = c->radius + c->pair->radius;
    double unpair_th = unpair_dist_coef*rad_sum;
    double distSq = 0;
    if (c->spr_nat_len < 2.0*c->radius) {
        c->spr_nat_len += DT_Cell*eps_L;
        c->pair->spr_nat_len = c->spr_nat_len;
        //fprintf(stdout, "nat len calced.\n");
        //fflush(stdout);
    }
    else if ((distSq = p_cell_dist_sq(c, c->pair))>unpair_th*unpair_th) {
        c->spr_nat_len = 0;
        c->pair->spr_nat_len = 0;
        c->pair->pair = nullptr;
        c->pair = nullptr;
        printf("unpaired. distSq:%lf\n", distSq);
    }

    //fprintf(stdout, "disperse end.\n");
}

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//    外部公開関数                                                            //
//                                                                            //
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
    //fprintf(stdout, "disperse parallel start.\n");
    //fflush(stdout);
    cman.other_foreach_parallel_native([](Cell*const RESTRICT c) {
        //auto&c = cman[i];
       
        
        if (c->pair != nullptr) {
            //fprintf(stdout, "pair exist.\n");
            //fflush(stdout);
            if (no_double_count(c, c->pair)) {
                //fprintf(stdout, "first pair confirmed.\n");
                //fflush(stdout);
                pair_disperse(c);
                //fprintf(stdout, "first pair proc end.\n");
                //fflush(stdout);
            }
            //fprintf(stdout, "pair pair proc end.\n");
            //fflush(stdout);
        }
    });
    //fprintf(stdout, "disperse parallel end.\n");
    //fflush(stdout);
}
/**
 *  強制的にカルシウム刺激を与える場合の状態の初期化を行う。
 *  @param cman CellManager
 *  @param [in] zzmax 細胞が存在する最大高さ(z座標)
 */
void initialize_sc(CellManager & cman, double zzmax)
{
    cman.other_foreach_parallel_native([zzmax](Cell*const RESTRICT c) {
        if (c->state == ALIVE&&zzmax - c->z() < 8 * cont::R_max) {
            c->agek = c->agek > cont::THRESH_DEAD ? c->agek : cont::THRESH_DEAD;
            c->state = AIR;
        }
    });
}
