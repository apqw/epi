#include "ext_stim.h"
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/blocked_range3d.h>
#include "cell.h"
#include "../../global.h"

inline real fB(real age, real B, bool cornif) {

   
    return (cornif&&age > pm->THRESH_DEAD - pm->DUR_ALIVE&&age <= pm->THRESH_DEAD + pm->DUR_DEAD ? 1 : 0) - pm->kb*B;
}

void calc_ext_stim(SwapData<Dyn3DArr<real>>& ext_stim, const Dyn3DArr<const Cell*>& cmap1, const Dyn3DArr<uint_fast8_t>& cmap2, real zzmax) {
   // static constexpr double DB = 0.0009;
   // using namespace cont;
    const int iz_bound = (int)((zzmax + pm->FAC_MAP*pm->R_max) *pm->inv_dz);
    auto& carr = ext_stim.first();
    auto& narr = ext_stim.second();
    const real dtc = pm->DT_Ca;
    const real db = pm->DB;
    const real idzSq = pm->inv_dz*pm->inv_dz;
    tbb::parallel_for(tbb::blocked_range3d<size_t>(0, pm->NX, 0, pm->NY, 0, iz_bound), [&](const tbb::blocked_range3d<size_t>& range) {

        for (size_t j = range.pages().begin(); j < range.pages().end(); ++j) {
            const int prev_x = per_x_prev_idx[j];
            const int next_x = per_x_next_idx[j];
            for (size_t k = range.rows().begin(); k < range.rows().end(); ++k) {
                const int prev_y = per_y_prev_idx[k];
                const int next_y = per_y_next_idx[k];
                for (size_t l = range.cols().begin(); l < range.cols().end(); ++l) {
                    const int prev_z = per_z_prev_idx[l], next_z = per_z_next_idx[l];
                    real dum_age = 0;
                    bool flg_cornified = false;

                    if (cmap1.at(j,k,l) != nullptr) {
                        if (cmap1.at(j, k, l)->state == DEAD || cmap1.at(j, k, l)->state == ALIVE) {
                            dum_age = cmap1.at(j, k, l)->agek;
                            flg_cornified = true;
                        }
                    }

                    auto& cext = carr.at(j, k, l);

                    narr.at(j, k, l) = cext
                        + dtc*(db *
                        (cmap2.at(prev_x,k,l) * (carr.at(prev_x,k,l) - cext)
                            + cmap2.at(j,prev_y,l) * (carr.at(j,prev_y,l) - cext)
                            + cmap2.at(j,k,prev_z) * (carr.at(j,k,prev_z) - cext)
                            + cmap2.at(j,k,next_z) * (carr.at(j,k,next_z) - cext)
                            + cmap2.at(j,next_y,l) * (carr.at(j,next_y,l) - cext)
                            + cmap2.at(next_x,k,l) * (carr.at(next_x,k,l) - cext)

                            ) *idzSq
                            + fB(dum_age, cext, flg_cornified));
                }
            }
        }
    });
    tbb::parallel_for(0, iz_bound, [&](int l) {
        for (size_t j = 0; j < pm->NX; j++) narr.at(j,pm->NY,l) = narr.at(j,0,l);
        for (size_t k = 0; k <= pm->NY; k++)narr.at(pm->NX,k,l) = narr.at(0,k,l);
    });

    ext_stim.swap();
}