#include "ca2p.h"
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/blocked_range3d.h>
#include "CellManager.h"
#include "../../global.h"
#include "../../misc/swapdata.h"
#include <cmath>
using GJSet = std::vector<std::vector<real>>;
inline void init_ca2p_avg(CellManager& cman) {
    cman.other_foreach([](Cell* const RESTRICT c) {
        auto& st = c->state;
        if (st == ALIVE || st == FIX || st == MUSUME) {
            c->ca2p_avg = 0;
        }
    });
}

////////////////////////////////////////////////////////////////////////////////
//                                                                          
//  Ca2+の反応項に関する定義
//
////////////////////////////////////////////////////////////////////////////////

inline real ca2p_by_ext_stim(real ext_stim) {

    return pm->kbc * ext_stim*ext_stim * pm->Cout / (pm->Hb + ext_stim*ext_stim);
}

inline real ca2p_into_storage(real ca2p_in_cell) {
    return pm->gamma * ca2p_in_cell / (pm->kg + ca2p_in_cell);
}

inline real ER_domain1_active(real ip3) {
    return (pm->mu0 + pm->mu1 * ip3 / (ip3 + pm->kmu));
}

inline real ER_domain2_active(real ca2p) {
    return (pm->para_b + pm->para_bb * ca2p / (pm->para_k1 + ca2p));
}

/**
*  Ca2+の反応項
*  @param [in] u Ca2+濃度
*  @param [in] ER_domain3_deactive 貯蔵庫から細胞質内への放出に対する不活性効果
*  @param [in] p IP3濃度
*  @param [in] B 細胞外刺激物質濃度
*/
inline real ca2p_reaction_factor(real u, real ER_domain3_deactive, real p, real B)
{

    real ca2p_from_storage = pm->k_flux * ER_domain1_active(p) * ER_domain2_active(u) * ER_domain3_deactive;

    return
        ca2p_from_storage
        - ca2p_into_storage(u)
        + pm->leak_from_storage
        + ca2p_by_ext_stim(B);

}


inline real calc_th(bool is_alive, real agek) {



    return is_alive ?
        pm->thgra + ((pm->thpri - pm->thgra)*real(0.5)) * (real(1.0) + tanh((pm->THRESH_SP - agek) / pm->delta_th)) :
        pm->thpri;
}

inline real calc_Kpa(bool is_alive, real agek) {

    return is_alive ?
        pm->Kgra + ((pm->Kpri - pm->Kgra)*real(0.5)) * (real(1.0) + tanh((pm->THRESH_SP - agek) / pm->delta_K)) :
        pm->Kpri;
}

inline real calc_IAG(bool is_alive, real agek) {
    return is_alive ?
        real(0.5)*(real(1.0) + tanh((agek - pm->THRESH_SP) / pm->delta_I)) :
        pm->iage_kitei;
}

inline real ex_inert_diff(real ca2p,real current_ex_inert, real _th) {
    return ((pm->para_k2*pm->para_k2 / (pm->para_k2*pm->para_k2 + ca2p*ca2p) - current_ex_inert) / _th);
}

inline real IP3_default_diff(real _Kpa, real a_avg, real current_IP3) {
    return (_Kpa*a_avg / (pm->H0 + a_avg) - pm->Kpp*current_IP3);
}

inline real fw(real diff, real w)
{
    return (real(1.) - w) + (real(-1.) + tanh((pm->wd - diff) / real(0.1))) / real(2.); //epsw0 == 0.1 //done
                                                            //-1. <-????
}

////////////////////////////////////////////////////////////////////////////////
//                                                                          
//  細胞内IP3,Ca2+の計算に関する定義
//
////////////////////////////////////////////////////////////////////////////////

/**
*  角層のIP3の値の次の値を計算する。
*  @attention CellManager側でswapしない限り現在の値は更新されない。
*/
inline void dead_IP3_calc(CellManager& cman,SwapData<std::vector<real>>& IP3Storage) {
    cman.other_foreach_parallel_native([&](Cell*const RESTRICT c) {
        if (c->state == DEAD) {

            int count = 0;
            real tmp = 0;
            c->connected_cell.foreach_with_index([&count, &c, &tmp,&IP3Storage](Cell* conn,size_t index) {
                if (conn->state == ALIVE) {
                    tmp += IP3Storage.first()[conn->get_index()];
                    count++;
                }
            });
            const size_t idx = c->get_index();
            tmp = pm->DT_Ca*(pm->dp*(tmp - count*IP3Storage.first()[idx]) - pm->Kpp*IP3Storage.first()[idx]);
            IP3Storage.second()[idx] = IP3Storage.first()[idx] + tmp;
        }
    });
}

/**
*  角層以外のIP3、カルシウム濃度の値の次の値を計算する。
*  @param cman CellManager
*  @param [in] ATP_first ATPの場の値
*  @param [in] ext_stim_first ext_stimの場の値
*  @attention CellManager側でswapしない限り現在の値は更新されない。
*/

inline real grid_avg8(const Dyn3DArr<real>& mat,int ix,int iy,int iz) {
    return real(0.125)*(mat.at(ix, iy, iz) +
        mat.at(ix + 1, iy, iz) + mat.at(ix, iy + 1, iz) + mat.at(ix, iy, iz + 1) +
        mat.at(ix, iy + 1, iz + 1) + mat.at(ix + 1, iy, iz + 1) + mat.at(ix + 1, iy + 1, iz) +
        mat.at(ix + 1, iy + 1, iz + 1));
}
inline void supra_calc(CellManager& cman, const Dyn3DArr<real>& ATP_first, const Dyn3DArr<real>& ext_stim_first,
    GJSet& GJStorage, SwapData<std::vector<real>>& IP3Storage, SwapData<std::vector<real>>& Ca2PStorage) {
    cman.other_foreach_parallel_native([&](Cell*const RESTRICT c) {
        auto& st = c->state;
        const size_t my_idx = c->get_index();
        if (st == ALIVE || st == FIX || st == MUSUME) {
            //aliasing
            int& ix = c->lat[0];
            int& iy = c->lat[1];
            int& iz = c->lat[2];

            const bool is_alive = c->state == ALIVE;
            const real _th = calc_th(is_alive, c->agek);
            const real _Kpa = calc_Kpa(is_alive, c->agek);
            const real IAGv = calc_IAG(is_alive, c->agek);
            const real myCa2P = Ca2PStorage.first()[my_idx];
            const real myIP3 = IP3Storage.first()[my_idx];

            real tmp_diffu = 0;
            real tmp_IP3 = 0;
            c->connected_cell.foreach_with_index([&](Cell*const RESTRICT conn,size_t conn_idx) {
                auto& st = conn->state;
                if (st == ALIVE || st == DEAD || st == FIX || st == MUSUME) {
                    const real conngj = GJStorage[my_idx][conn_idx];
                    const real connCa2P = Ca2PStorage.first()[conn->get_index()];
                    
                    const real connIP3 = IP3Storage.first()[conn->get_index()];
                    
                    tmp_diffu += conngj * (connCa2P - myCa2P);
                    tmp_IP3 += conngj*(connIP3 - myIP3);
                    if (st == ALIVE) {
                        GJStorage[my_idx][conn_idx] += pm->DT_Ca*fw(fabs(connCa2P - myCa2P), conngj);
                    }
                }
            });

            c->diff_u = ca2p_reaction_factor(myCa2P, c->ex_inert, myIP3, grid_avg8(ext_stim_first, ix, iy, iz)) + pm->ca2p_du*IAGv*tmp_diffu;

            IP3Storage.second()[my_idx] = myIP3 + pm->DT_Ca*(IP3_default_diff(_Kpa, grid_avg8(ATP_first, ix, iy, iz), myIP3) + pm->dp*IAGv*tmp_IP3);
            const real cs = myCa2P + pm->DT_Ca*c->diff_u;
            Ca2PStorage.second()[my_idx] = cs;
            c->ca2p_avg += cs;

            c->ex_inert += pm->DT_Ca*ex_inert_diff(myCa2P, c->ex_inert, _th); //update last
        }

    });
}

////////////////////////////////////////////////////////////////////////////////
//                                                                          
//  場に関するマップ作成
//
////////////////////////////////////////////////////////////////////////////////

/**
*  場全体に渡るマップを作成。平均化ごとに1度だけ計算する。
*  @param [out] air_stim_flg 各座標について、空気による刺激が有効かどうかのマップ
*  @param [out] cell_diffu_map 各座標について、その位置にある細胞の、カルシウム差分への参照を割り当てる。
*  @param [in] cmap1 各座標について、細胞が存在するかどうかのマップ
*  @param [in] iz_bound 細胞が存在する上限(z座標) (半径も含む)
*  @param [in] default_diffu_ptr 細胞が存在しない座標におけるcell_diffu_mapの参照先。
*/
void init_ca2p_map(Dyn3DArr<uint_fast8_t>& air_stim_flg, Dyn3DArr<const real*>& cell_diffu_map, const Dyn3DArr<const Cell*>& cmap1, int iz_bound, const real*const default_diffu_ptr) {
    tbb::parallel_for(tbb::blocked_range3d<int>(0, pm->NX, 0, pm->NY, 0, iz_bound), [&](const tbb::blocked_range3d<int>& range) {
        for (int j = range.pages().begin(); j != range.pages().end(); ++j) {
            for (int k = range.rows().begin(); k < range.rows().end(); k++) {
                for (int l = range.cols().begin(); l < range.cols().end(); l++) {
                    const real* tmp = default_diffu_ptr;
                    uint_fast8_t asf = 0;
                    auto& c = cmap1.at(j, k, l);
                    if (c != nullptr) {

                        /*
                        ptr to diffu
                        */
                        auto& st = c->state;
                        if (st == ALIVE || st == FIX || st == MUSUME) {
                            tmp = &(c->diff_u);
                        }


                        /*
                        connected with air
                        */
                        size_t c_count = c->connected_cell.size();
                        for (size_t cc = 0; cc < c_count; ++cc) {
                            if (c->connected_cell[cc]->state == AIR) {
                                asf = 1;
                                break;
                            }
                        }
                    }

                    cell_diffu_map.at(j,k,l) = tmp;
                    air_stim_flg.at(j,k,l) = asf;
                }
            }
        }
    });

}


////////////////////////////////////////////////////////////////////////////////
//                                                                          
//  ATPの計算に関する定義
//
////////////////////////////////////////////////////////////////////////////////

inline real fa(real diffu, real A) {
    return pm->STIM11*min0(diffu) - A*pm->Kaa;
}

/**
* ATPの次の値を計算。
* @param ATP ATPの場の値(現在の値と次の値の格納先のセット)
* @param [in] cell_diffu_map 各座標に存在する細胞のカルシウム差分への参照
* @param [in] air_stim_flg 各座標について、空気による刺激が有効かどうかのマップ
* @param [in] cmap2 差分用マップ
* @param [in] iz_bound 細胞が存在する上限(z座標) (半径も含む)
* @attention ATPをswapしない限り現在の値は更新されない。
*/
inline void ATP_refresh(SwapData<Dyn3DArr<real>>& ATP, const Dyn3DArr<const real*>& cell_diffu_map, const Dyn3DArr<uint_fast8_t>& air_stim_flg, const Dyn3DArr<uint_fast8_t>& cmap2, int iz_bound) {

    auto& carr = ATP.first();
    auto&& narr = ATP.second();
    auto& cell_map2 = cmap2;
    tbb::parallel_for(tbb::blocked_range3d<int>(0, pm->NX, 0, pm->NY, 0, iz_bound), [&](const tbb::blocked_range3d<int>& range) {
        for (int j = range.pages().begin(); j != range.pages().end(); ++j) {
            const int prev_x = precalc_per_prev_x()[j];
            const int next_x = precalc_per_next_x()[j];
            for (int k = range.rows().begin(); k != range.rows().end(); ++k) {
                const int prev_y = precalc_per_prev_y()[k];
                const int next_y = precalc_per_next_y()[k];
                for (int l = range.cols().begin(); l != range.cols().end(); ++l) {
                    const int prev_z = precalc_per_prev_z()[l];
                    const int next_z = precalc_per_next_z()[l];
                    real& catp = carr.at(j,k,l);
                    narr.at(j,k,l) = catp + pm->DT_Ca*(pm->Da * (cell_map2.at(prev_x,k,l) * (carr.at(prev_x,k,l) - catp)
                        + cell_map2.at(next_x,k,l) * (carr.at(next_x,k,l) - catp)
                        + cell_map2.at(j,prev_y,l) * (carr.at(j,prev_y,l) - catp)
                        + cell_map2.at(j,next_y,l) * (carr.at(j,next_y,l) - catp)
                        + cell_map2.at(j,k,prev_z) * (carr.at(j,k,prev_z) - catp)
                        + cell_map2.at(j,k,next_z) * (carr.at(j,k,next_z) - catp)) *pm->inv_dx*pm->inv_dx
                        + fa(*cell_diffu_map.at(j,k,l), catp) + air_stim_flg.at(j,k,l) * pm->AIR_STIM);
                }
            }
        }
    });

    for (int l = 0; l <= iz_bound; l++) {
        for (size_t j = 0; j <= pm->NX; j++) narr.at(j, pm->NY,l) = narr.at(j,0,l);
        for (size_t k = 0; k <= pm->NY; k++) narr.at(pm->NX,k,l) = narr.at(0,k,l);
    }
}
/**
*  値の更新。
*  swap等を行う。
*/
inline void update_values(CellManager& cman, SwapData<Dyn3DArr<real>>& ATP, SwapData<std::vector<real>>& IP3Storage, SwapData<std::vector<real>>& Ca2PStorage) {
    Ca2PStorage.swap();
    IP3Storage.swap();
    cman.other_foreach([](Cell* const RESTRICT c) {
        c->diff_u = 0;
    });
    ATP.swap();
}

/**
*  細胞にカルシウム濃度の平均値をセットする。
*/
inline void set_cell_ca2p(CellManager& cman, const std::vector<real>& Ca2PStorage_current, const std::vector<real>& IP3Storage_current) {
    cman.other_foreach_parallel_native([&](Cell*const RESTRICT c) {
        auto& st = c->state;
        c->IP3 = IP3Storage_current[c->get_index()];
        if (st == ALIVE || st == FIX || st == MUSUME) {
            c->ca2p_avg /= pm->Ca_ITR;
            c->ca2p = Ca2PStorage_current[c->get_index()];
        }
        else {
            c->ca2p_avg = 0;
            c->ca2p = 0;
        }
    });
}
/**
*  細胞のカルシウム濃度を計算。
* @param cman CellManager
* @param ATP ATPの場の値(現在の値と次の値の格納先のセット)
* @param [in] ext_stim_first ext_stimの場の値
* @param [in] cmap1 細胞が存在するかどうかのマップ
* @param [in] cmap2 差分用マップ
* @param [in] zzmax 細胞が存在する上限(z座標) (半径は含まない)
*/

void calc_ca2p(CellManager& cman, SwapData<Dyn3DArr<real>>& ATP, const Dyn3DArr<real>& ext_stim_first, const Dyn3DArr<const Cell*>& cmap1, const Dyn3DArr<uint_fast8_t>& cmap2, real zzmax) {
    //using namespace cont;
    init_ca2p_avg(cman);


   const int iz_bound = (int)((zzmax + pm->FAC_MAP*pm->R_max) *pm->inv_dz);
    static Dyn3DArr<uint_fast8_t> air_stim_flg(pm->NX+1,pm->NY+1,pm->NZ+1);
    static Dyn3DArr<const real*> cell_diffu_map(pm->NX + 1, pm->NY + 1, pm->NZ + 1);
    SwapData<std::vector<real>> IP3Storage{ cman.size() };
    SwapData<std::vector<real>> Ca2PStorage{ cman.size() };
    GJSet GJStorage{ cman.size() };


    cman.all_foreach_parallel_native([&](const Cell* c) {
        const size_t idx = c->get_index();
        IP3Storage.first()[idx] = c->IP3;
        Ca2PStorage.first()[idx] = c->ca2p;
        GJStorage[idx].resize(c->connected_cell.size(),pm->gj_init);
    });

    const real dummy_diffu = 0;
    init_ca2p_map(air_stim_flg, cell_diffu_map, cmap1, iz_bound, &dummy_diffu);

    for (size_t cstp = 0; cstp < pm->Ca_ITR; cstp++) {

        dead_IP3_calc(cman,IP3Storage);
        supra_calc(cman, ATP.first(), ext_stim_first,GJStorage,IP3Storage,Ca2PStorage);
        ATP_refresh(ATP, cell_diffu_map, air_stim_flg, cmap2, iz_bound);
        update_values(cman, ATP,IP3Storage,Ca2PStorage);

    }

    set_cell_ca2p(cman,Ca2PStorage.first(),IP3Storage.first());
}

