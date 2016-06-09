#include "ca2p.h"
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/blocked_range3d.h>
#include "utils.h"
#include "cell.h"
#include "cellmanager.h"
#include <cmath>

template<typename ArrTy>
inline auto grid_avg8(const ArrTy& grd,int ix,int iy,int iz) {
	return 0.125*(grd[ix][iy][iz] + grd[ix + 1][iy][iz] + grd[ix][iy + 1][iz]
		+ grd[ix][iy][iz + 1] + grd[ix + 1][iy + 1][iz] + grd[ix + 1][iy][iz + 1]
		+ grd[ix][iy + 1][iz + 1] + grd[ix + 1][iy + 1][iz + 1]);
}
inline void init_ca2p_avg(CellManager& cman) {
    cman.other_foreach([](Cell* const RESTRICT c) {
		auto& st = c->state;
		if (st == ALIVE || st == FIX || st == MUSUME) {
			c->ca2p_avg = 0;
		}
	});
}

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
			c->IP3 = c->IP3() + tmp;
		}
	});
}
inline double fu(double u, double v, double p, double B)
{
	//  static double result_debug=0.0;
	using namespace cont;
	constexpr double c_gamma = cont::gamma;
	return kf * (mu0 + mu1 * p / (p + kmu)) * (para_b + para_bb * u / (para_k1 + u)) * v
		- c_gamma * u / (kg + u) + beta + kbc * B*B * Cout / (Hb + B*B);

}

inline double ALIVE_th(double agek) {
	using namespace cont;
	return thgra + ((thpri - thgra)*0.5) * (1.0 + tanh((THRESH_SP - agek) / delta_th));
}

inline double ALIVE_Kpa(double agek) {
	using namespace cont;
	return Kgra + ((Kpri - Kgra)*0.5) * (1.0 + tanh((THRESH_SP - agek) / delta_K));
}

inline double ALIVE_IAG(double agek) {
	using namespace cont;
	return 0.5*(1.0 + tanh((agek - THRESH_SP) / delta_I));
}

inline double ex_inert_diff(double ca2p, double current_ex_inert,double _th) {
	using namespace cont;
	return ((para_k2*para_k2 / (para_k2*para_k2 + ca2p*ca2p) - current_ex_inert) / _th);
}

inline double IP3_default_diff(double _Kpa,double a_avg,double current_IP3) {
	using namespace cont;
	return (_Kpa*a_avg / (H0 + a_avg) - Kpp*current_IP3);
}

inline double fw(double diff, double w)
{
	using namespace cont;
	return (1. - w) + (-1. + tanh((wd - diff) / 0.1)) / 2.; //epsw0 == 0.1 //done
															//-1. <-????
}

inline void supra_calc(CellManager& cman,const FArr3D<double>& ATP_first, const  FArr3D<double>& ext_stim_first) {
	using namespace cont;
	cman.other_foreach_parallel_native([&](Cell*const RESTRICT c){
		auto& st = c->state;
		if (st == ALIVE || st == FIX || st == MUSUME) {
			//aliasing
			int& ix = c->lat[0];
			int& iy = c->lat[1];
			int& iz = c->lat[2];

			const bool is_alive = c->state == ALIVE;
			const double _th = is_alive? ALIVE_th(c->agek) :thpri;
			const double _Kpa = is_alive? ALIVE_Kpa(c->agek): Kpri;
			const double IAGv = is_alive? ALIVE_IAG(c->agek):iage_kitei;

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

			c->diff_u= fu(c->ca2p(), c->ex_inert, c->IP3(), grid_avg8(ext_stim_first(), ix, iy, iz))+ca2p_du*IAGv*tmp_diffu;
			c->IP3 = c->IP3() + DT_Ca*(IP3_default_diff(_Kpa, grid_avg8(ATP_first(), ix, iy, iz), c->IP3()) + dp*IAGv*tmp_IP3);
			const double cs = c->ca2p() + DT_Ca*c->diff_u;
			c->ca2p = cs;
			c->ca2p_avg += cs;

			c->ex_inert += DT_Ca*ex_inert_diff(c->ca2p(), c->ex_inert, _th); //update last
		}
	
	});
}


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
						int c_count = c->connected_cell.size();
						for (int cc = 0; cc < c_count; ++cc) {
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
inline double fa(double diffu, double A) {
	using namespace cont;
	return STIM11*min0(diffu) - A*Kaa;
}
inline void ATP_refresh(SwapData<FArr3D<double>>& ATP, const RawArr3D<const double*>& cell_diffu_map, const RawArr3D<uint_fast8_t>& air_stim_flg, const FArr3D<uint_fast8_t>& cmap2,int iz_bound) {
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

inline void update_values(CellManager& cman, SwapData<FArr3D<double>>& ATP) {
	ca2p_swap(cman);
	IP3_swap(cman);
    cman.other_foreach([](Cell* const RESTRICT c) {
		c->diff_u = 0;
	});
	ATP.swap();
}

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

void calc_ca2p(CellManager& cman, SwapData<FArr3D<double>>& ATP,const FArr3D<double>& ext_stim_first,const FArr3D<const Cell*>& cmap1,const FArr3D<uint_fast8_t>& cmap2, double zzmax)
{
	using namespace cont;
	init_ca2p_avg(cman);


	int iz_bound = (int)((zzmax + FAC_MAP*R_max) / dz);
	static RawArr3D<uint_fast8_t> air_stim_flg;
	static RawArr3D<const double*> cell_diffu_map;
	const double dummy_diffu = 0;
	init_ca2p_map(air_stim_flg, cell_diffu_map, cmap1, iz_bound, &dummy_diffu);

	for (int cstp = 0; cstp < Ca_ITR; cstp++) {
		
		dead_IP3_calc(cman);
		supra_calc(cman, ATP.first(), ext_stim_first);
		ATP_refresh(ATP, cell_diffu_map, air_stim_flg, cmap2, iz_bound);
		update_values(cman, ATP);

	}

	set_cell_ca2p(cman);
}
