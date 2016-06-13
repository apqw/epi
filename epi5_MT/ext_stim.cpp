#include "ext_stim.h"
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/blocked_range3d.h>
#include "utils.h"
#include "cell.h"

inline double fB(double age, double B, bool cornif) {
	using namespace cont;
	return (cornif&&age > THRESH_DEAD - DUR_ALIVE&&age <= THRESH_DEAD + DUR_DEAD ? 1 : 0) - kb*B;
}

void calc_ext_stim(SwapData<FArr3D<double>>& ext_stim,const FArr3D<const Cell*>& cmap1, const FArr3D<cmask_ty>& cmap2, double zzmax)
{
	using namespace cont;
	const int iz_bound = (int)((zzmax + FAC_MAP*R_max) *inv_dz);
	static int a_prev_z[cont::NZ];
	static int a_next_z[cont::NZ];
	for (int l = 0; l < iz_bound; l++) {
		int prev_z = 0, next_z = 0;
		if (l == 0) {
			prev_z = 1;
		}
		else {
			prev_z = l - 1;
		}
		if (l == NZ) {
			next_z = NZ - 1;
		}
		else {
			next_z = l + 1;
		}
		a_prev_z[l] = prev_z;
		a_next_z[l] = next_z;
	}
    auto&& carr = ext_stim.first()();
    auto&& narr = ext_stim.second()();

	tbb::parallel_for(tbb::blocked_range3d<size_t>(0, NX, 0, NY, 0, iz_bound), [&](const tbb::blocked_range3d<size_t>& range) {
		
		for (size_t j = range.pages().begin(); j < range.pages().end(); ++j) {
			const int prev_x = per_x_prev_idx[j];
			const int next_x = per_x_next_idx[j];
			for (size_t k = range.rows().begin(); k < range.rows().end(); ++k) {
				const int prev_y = per_y_prev_idx[k];
				const int next_y = per_y_next_idx[k];
                for (size_t l = range.cols().begin(); l < range.cols().end(); ++l) {
					const int prev_z = a_prev_z[l], next_z = a_next_z[l];
					double dum_age = 0;
					bool flg_cornified = false;

					if (cmap1()[j][k][l] != nullptr) {
						if (cmap1()[j][k][l]->state_mask()&(DEAD_M|ALIVE_M)) {
							dum_age = cmap1()[j][k][l]->agek;
							flg_cornified = true;
						}
					}
					
					auto& cext = carr[j][k][l];

					narr[j][k][l] = cext
						+ DT_Ca*(DB *
							(cmap2()[prev_x][k][l] * (carr[prev_x][k][l] - cext)
								+ cmap2()[j][prev_y][l] * (carr[j][prev_y][l] - cext)
								+ cmap2()[j][k][prev_z] * (carr[j][k][prev_z] - cext)
								+ cmap2()[j][k][next_z] * (carr[j][k][next_z] - cext)
								+ cmap2()[j][next_y][l] * (carr[j][next_y][l] - cext)
								+ cmap2()[next_x][k][l] * (carr[next_x][k][l] - cext)
								
								) *inv_dz*inv_dz
							+ fB(dum_age, cext, flg_cornified));
				}
			}
		}
	});

    for (int l = 0; l <= iz_bound; l++) {
        for (size_t j = 0; j < NX; j++) narr[j][NY][l] = narr[j][0][l];
        for (size_t k = 0; k <= NY; k++)narr[NX][k][l] = narr[0][k][l];
	}

	ext_stim.swap();
}

