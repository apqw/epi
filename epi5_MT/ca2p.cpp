#include "ca2p.h"
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/blocked_range3d.h>
#include "utils.h"
#include "cell.h"
#include "cellmanager.h"
#include <cmath>
/**
 *  @file カルシウム濃度、ATPの計算
 */

/** カルシウムの平均化に使うイテレーション回数 */
static constexpr unsigned Ca_ITR = (int)(cont::Ca_avg_time / cont::DT_Ca);

/** @todo 分かりやすい命名 */                                      
static constexpr double Kpp		= 0.3;
/** @todo 分かりやすい命名 */
static constexpr double dp		= 0.1;


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

inline double ca2p_by_ext_stim(double ext_stim) {
    static constexpr double kbc = 0.4*1.2;
    static constexpr double Hb = 0.01;
    static constexpr double Cout = 1.0;

    return kbc * ext_stim*ext_stim * Cout / (Hb + ext_stim*ext_stim);
}

inline double ca2p_into_storage(double ca2p_in_cell) {
    static constexpr double kg = 0.1;
    static constexpr double gamma = 2.0;
    return gamma * ca2p_in_cell / (kg + ca2p_in_cell);
}

inline double ER_domain1_active(double ip3) {
    static constexpr double mu0 = 0.567;
    static constexpr double mu1 = 0.1;
    static constexpr double kmu = 0.05;
    return (mu0 + mu1 * ip3 / (ip3 + kmu));
}

inline double ER_domain2_active(double ca2p) {
    static constexpr double para_b = 0.11;
    static constexpr double para_bb = 0.89;
    static constexpr double para_k1 = 0.7;
    return (para_b + para_bb * ca2p / (para_k1 + ca2p));
}

/** 
 *  Ca2+の反応項
 *  @param [in] u Ca2+濃度
 *  @param [in] ER_domain3_deactive 貯蔵庫から細胞質内への放出に対する不活性効果
 *  @param [in] p IP3濃度
 *  @param [in] B 細胞外刺激物質濃度
 */
inline double ca2p_reaction_factor(double u, double ER_domain3_deactive, double p, double B)
{
	static constexpr double k_flux		= 8.1;

	static constexpr double beta_zero	= 0.02;
	static constexpr double CA_OUT		= 1.0;
	static constexpr double leak_from_storage		= CA_OUT*beta_zero;

    double ca2p_from_storage = k_flux * ER_domain1_active(p) * ER_domain2_active(u) * ER_domain3_deactive;
	
	return 
          ca2p_from_storage
		- ca2p_into_storage(u) 
        + leak_from_storage 
        + ca2p_by_ext_stim(B);

}

/** @todo 分かりやすい命名 */
inline double calc_th(bool is_alive,double agek) {
	
	static constexpr double thgra		= 0.2;
	static constexpr double delta_th	= 1.0;
	static constexpr double thpri = 1.0;
	using namespace cont;
	
	
	return is_alive ?
		thgra + ((thpri - thgra)*0.5) * (1.0 + tanh((THRESH_SP - agek) / delta_th)) :
		thpri;
}

/** @todo 分かりやすい命名 */
inline double calc_Kpa(bool is_alive,double agek) {
	
	static constexpr double kpa			= 4.0;
	static constexpr double Kgra		= kpa;
	static constexpr double delta_K		= 1.0;
	static constexpr double Kpri = 6.0;

	using namespace cont;
	return is_alive ?
		Kgra + ((Kpri - Kgra)*0.5) * (1.0 + tanh((THRESH_SP - agek) / delta_K)) :
		Kpri;
}

/** @todo 分かりやすい命名 */
inline double calc_IAG(bool is_alive,double agek) {
	static constexpr double delta_I = 1.5;
	static constexpr double iage_kitei = 0;//ok
	using namespace cont;
	return is_alive ?
		0.5*(1.0 + tanh((agek - THRESH_SP) / delta_I)) :
		iage_kitei;
}

/** @todo 分かりやすい命名 */
inline double ex_inert_diff(double ca2p, double current_ex_inert,double _th) {
	static constexpr double para_k2 = 0.7;
	using namespace cont;
	return ((para_k2*para_k2 / (para_k2*para_k2 + ca2p*ca2p) - current_ex_inert) / _th);
}

/** @todo 分かりやすい命名 */
inline double IP3_default_diff(double _Kpa,double a_avg,double current_IP3) {
	static constexpr double H0 = 0.5;

	using namespace cont;
	return (_Kpa*a_avg / (H0 + a_avg) - Kpp*current_IP3);
}

/** @todo 分かりやすい命名 */
inline double fw(double diff, double w)
{
	static constexpr double wd = 0.1;//ok
	using namespace cont;
	return (1. - w) + (-1. + tanh((wd - diff) / 0.1)) / 2.; //epsw0 == 0.1 //done
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
inline void dead_IP3_calc(CellManager& cman) {
    using namespace cont;
    cman.other_foreach_parallel_native([&](Cell*const RESTRICT c) {
        if (c->state == DEAD) {

            int count = 0;
            double tmp = 0;
            c->connected_cell.foreach([&count, &c, &tmp](Cell* conn) {
                if (conn->state == ALIVE) {
                    tmp += conn->IP3();
                    count++;
                }
            });
            tmp = DT_Ca*(dp*(tmp - count*c->IP3()) - Kpp*c->IP3());
            c->IP3.set_next(c->IP3() + tmp);
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
inline void supra_calc(CellManager& cman,const FArr3D<double>& ATP_first, const  FArr3D<double>& ext_stim_first) {
	static constexpr double ca2p_du = 0.01;
	
	using namespace cont;
	cman.other_foreach_parallel_native([&](Cell*const RESTRICT c){
		auto& st = c->state;
		if (st == ALIVE || st == FIX || st == MUSUME) {
			//aliasing
			int& ix = c->lat[0];
			int& iy = c->lat[1];
			int& iz = c->lat[2];

			const bool is_alive = c->state == ALIVE;
			const double _th = calc_th(is_alive,c->agek);
			const double _Kpa = calc_Kpa(is_alive,c->agek);
			const double IAGv = calc_IAG(is_alive, c->agek);

			double tmp_diffu = 0;
			double tmp_IP3 = 0;
            c->connected_cell.foreach([&c, &tmp_diffu, &tmp_IP3, IAGv](Cell*const RESTRICT conn) {
				auto& st = conn->state;
				if (st==ALIVE||st==DEAD||st==FIX||st==MUSUME) {
					tmp_diffu += c->gj.at(conn)()*(conn->ca2p() - c->ca2p());
					tmp_IP3 += c->gj.at(conn)()*(conn->IP3() - c->IP3());
					if (st == ALIVE) {
						c->gj.at(conn)() += DT_Ca*fw(fabs(conn->ca2p() - c->ca2p()), c->gj.at(conn)());
					}
				}
			});

			c->diff_u= ca2p_reaction_factor(c->ca2p(), c->ex_inert, c->IP3(), grid_avg8(ext_stim_first(), ix, iy, iz))+ca2p_du*IAGv*tmp_diffu;
			c->IP3.set_next(c->IP3() + DT_Ca*(IP3_default_diff(_Kpa, grid_avg8(ATP_first(), ix, iy, iz), c->IP3()) + dp*IAGv*tmp_IP3));
			const double cs = c->ca2p() + DT_Ca*c->diff_u;
			c->ca2p.set_next(cs);
			c->ca2p_avg += cs;

			c->ex_inert += DT_Ca*ex_inert_diff(c->ca2p(), c->ex_inert, _th); //update last
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
void init_ca2p_map(RawArr3D<uint_fast8_t>& air_stim_flg, RawArr3D<const double*>& cell_diffu_map,const FArr3D<const Cell*>& cmap1,int iz_bound,const double*const default_diffu_ptr) {
	using namespace cont;

	tbb::parallel_for(tbb::blocked_range3d<int>(0, NX, 0, NY, 0, iz_bound), [&](const tbb::blocked_range3d<int>& range) {
		for (int j = range.pages().begin(); j != range.pages().end(); ++j) {
			for (int k = range.rows().begin(); k < range.rows().end(); k++) {
				for (int l = range.cols().begin(); l < range.cols().end(); l++) {
					const double* tmp = default_diffu_ptr;
					uint_fast8_t asf = 0;
					auto& c = cmap1()[j][k][l];
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

					cell_diffu_map[j][k][l] = tmp;
					air_stim_flg[j][k][l] = asf;
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

/** @todo 分かりやすい命名 */
inline double fa(double diffu, double A) {
	static constexpr double STIM11 = 0.002;//ok
	static constexpr double Kaa = 0.5;//ok
	using namespace cont;
	return STIM11*min0(diffu) - A*Kaa;
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
inline void ATP_refresh(SwapData<FArr3D<double>>& ATP, const RawArr3D<const double*>& cell_diffu_map, const RawArr3D<uint_fast8_t>& air_stim_flg, const FArr3D<cmask_ty>& cmap2,int iz_bound) {
	static constexpr double Da = 1.0;
	static constexpr double AIR_STIM = 0.1;


	using namespace cont;
	auto&& carr = ATP.first()();
	auto&& narr = ATP.second()();
	auto& cell_map2 = cmap2();
	tbb::parallel_for(tbb::blocked_range3d<int>(0, NX, 0, NY, 0, iz_bound), [&](const tbb::blocked_range3d<int>& range) {
		for (int j = range.pages().begin(); j != range.pages().end(); ++j) {
			const int prev_x = precalc_per_prev_x()[j];
			const int next_x = precalc_per_next_x()[j];
			for (int k = range.rows().begin(); k != range.rows().end(); ++k) {
				const int prev_y = precalc_per_prev_y()[k];
				const int next_y = precalc_per_next_y()[k];
				for (int l = range.cols().begin(); l != range.cols().end(); ++l) {
					const int prev_z = precalc_per_prev_z()[l];
					const int next_z = precalc_per_next_z()[l];
					double& catp = carr[j][k][l];
					narr[j][k][l] = catp + DT_Ca*(Da * (cell_map2[prev_x][k][l] * (carr[prev_x][k][l] - catp)
						+ cell_map2[next_x][k][l] * (carr[next_x][k][l] - catp)
						+ cell_map2[j][prev_y][l] * (carr[j][prev_y][l] - catp)
						+ cell_map2[j][next_y][l] * (carr[j][next_y][l] - catp)
						+ cell_map2[j][k][prev_z] * (carr[j][k][prev_z] - catp)
						+ cell_map2[j][k][next_z] * (carr[j][k][next_z] - catp)) *inv_dx*inv_dx
						+ fa(*cell_diffu_map[j][k][l], catp) + air_stim_flg[j][k][l] * AIR_STIM);
				}
			}
		}
	});

    for (int l = 0; l <= iz_bound; l++) {
        for (size_t j = 0; j <= NX; j++) narr[j][NY][l] = narr[j][0][l];
        for (size_t k = 0; k <= NY; k++) narr[NX][k][l] = narr[0][k][l];
	}
}
/**
 *  値の更新。
 *  swap等を行う。
 */
inline void update_values(CellManager& cman, SwapData<FArr3D<double>>& ATP) {
	ca2p_swap(cman);
	IP3_swap(cman);
    cman.other_foreach([](Cell* const RESTRICT c) {
		c->diff_u = 0;
	});
	ATP.swap();
}

/**
 *  細胞にカルシウム濃度の平均値をセットする。
 */
inline void set_cell_ca2p(CellManager& cman) {
	using namespace cont;
    cman.other_foreach_parallel_native([](Cell*const RESTRICT c) {
		auto& st = c->state;
		if (st==ALIVE||st==FIX||st==MUSUME) {
			c->ca2p_avg /= Ca_ITR;
		}
		else {
			c->ca2p_avg=0;
			c->ca2p._set(0);
		}
	});
}

///////////////////////////////////////////////////////////////

/**
 *  細胞のカルシウム濃度を計算。
 * @param cman CellManager
 * @param ATP ATPの場の値(現在の値と次の値の格納先のセット)
 * @param [in] ext_stim_first ext_stimの場の値
 * @param [in] cmap1 細胞が存在するかどうかのマップ
 * @param [in] cmap2 差分用マップ
 * @param [in] zzmax 細胞が存在する上限(z座標) (半径は含まない)
 */
void calc_ca2p(CellManager& cman, SwapData<FArr3D<double>>& ATP,const FArr3D<double>& ext_stim_first,const FArr3D<const Cell*>& cmap1,const FArr3D<cmask_ty>& cmap2, double zzmax)
{
	using namespace cont;
	init_ca2p_avg(cman);


	int iz_bound = (int)((zzmax + FAC_MAP*R_max) / dz);
	static RawArr3D<uint_fast8_t> air_stim_flg;
	static RawArr3D<const double*> cell_diffu_map;
	const double dummy_diffu = 0;
	init_ca2p_map(air_stim_flg, cell_diffu_map, cmap1, iz_bound, &dummy_diffu);

    for (size_t cstp = 0; cstp < Ca_ITR; cstp++) {
		
		dead_IP3_calc(cman);
		supra_calc(cman, ATP.first(), ext_stim_first);
		ATP_refresh(ATP, cell_diffu_map, air_stim_flg, cmap2, iz_bound);
		update_values(cman, ATP);

	}

	set_cell_ca2p(cman);
}
