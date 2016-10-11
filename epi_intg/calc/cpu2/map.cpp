
#include "map.h"
#include "../../utils.h"
#include "../../global.h"
#include "CellManager.h"

#include <vector>
#include <cstdint>
/**
 *  細胞の格子点としての座標を計算
 */
inline void set_lattice(Cell*const RESTRICT c) {
	c->lat[0] = precalc_lat_x()[(int)(c->x()*pm->inv_dx) + pm->NX];
	c->lat[1] = precalc_lat_y()[(int)(c->y()*pm->inv_dy) + pm->NY];
	c->lat[2] = precalc_lat_z()[(int)(c->z()*pm->inv_dz) + pm->NZ];
}

void setup_map_lat(CellManager & cman, Dyn3DArr<const Cell*>& cmap1, Dyn3DArr<uint_fast8_t>& cmap2)
{

	const real Rad = pm->R_max;

	const int irx = (int)(Rad*pm->inv_dx);
	const int iry = (int)(Rad*pm->inv_dy);
	const int irz = (int)(Rad*pm->inv_dz);


	tbb::parallel_for(tbb::blocked_range<int>(0, cmap1.NX), [&cmap1,&cmap2](const tbb::blocked_range< int >& range) {
		cmap1.x_range_memset(range.begin(),range.size(),0);
		cmap2.x_range_memset(range.begin(),range.size(),0);
	});
	cman.all_foreach_parallel_native([](Cell*const RESTRICT c) {
		set_lattice(c);
	});
	cman.all_foreach_parallel_native([&](const Cell*const c) {

		auto& clat = c->lat;
		const bool afm_cell = c->state == ALIVE || c->state == FIX || c->state == MUSUME;

		 if (!afm_cell) {
			const int _irx = irx;
			const int _iry = iry;
			const int _irz = irz;
			const real normal_radSq = c->radius*c->radius;



			static thread_local std::vector<real> a_diffySq_n(_iry * 2 + 1);
			static thread_local std::vector<int> a_ipy_n(_iry * 2 + 1);
            const int ymax = clat[1] + _iry; const int ymin = clat[1] - _iry;
			for (int l = ymin, yc = 0; l <= ymax; l++) {
				const real my = l * pm->dy;
				a_ipy_n[yc] = precalc_lat_y()[l + pm->NY];
				const real diffy = p_diff_y(my, c->y());
				a_diffySq_n[yc] = diffy*diffy;
				yc++;
			}

            static thread_local std::vector<real> a_diffzSq_n(_irz * 2 + 1);
            const int z_b = clat[2] - _irz;
            const int z_u = clat[2] + _irz;
            const int zmin = z_b >= 1 ? z_b : 1;
            const int zmax = z_u < (int)pm->NZ ? z_u : (int)pm->NZ - 1;
			for (int m = zmin, zc = 0; m <= zmax; m++) {
				const real mz = m * pm->dz;
				const real diffz = mz - c->z();
				a_diffzSq_n[zc] = diffz*diffz;
				zc++;
			}
			//yc = 0; zc = 0;
            const int xmax = clat[0] + _irx; const int xmin = clat[0] - _irx;
			for (int k = xmin; k <= xmax; k++) {
				const real mx = k * pm->dx;
				const int ipx = precalc_lat_x()[k + pm->NX];
				const real diffx = p_diff_x(mx, c->x());
				const real diffxSq = diffx*diffx;
				for (int yc = 0; yc < _iry * 2 + 1; yc++) {
					const int ipy = a_ipy_n[yc];
					const real tmps = diffxSq + a_diffySq_n[yc];
					for (int zc = 0, ipz = zmin; ipz <= zmax; ipz++, zc++) {
						const real distSq = tmps + a_diffzSq_n[zc];
						if (distSq < normal_radSq) {
								cmap1.at(ipx,ipy,ipz) = c;
						}
					}
				}
			}
		}else {
			 const int _irx = 2 * irx;
			 const int _iry = 2 * iry;
			 const int _irz = 2 * irz;
			 //double diffx, diffxSq, diffy, diffz, distSq;
			 const real crad = pm->FAC_MAP*c->radius;
			 const real cradSq = crad*crad;
			 const real normal_radSq = c->radius*c->radius;
			 const int xmax = clat[0] + _irx; const int xmin = clat[0] - _irx;
			 const int ymax = clat[1] + _iry; const int ymin = clat[1] - _iry;
			 const int z_b = clat[2] - _irz;
			 const int z_u = clat[2] + _irz;
			 const int zmin = z_b >= 1 ? z_b : 1;
			 const int zmax = z_u < (int)pm->NZ ? z_u : (int)pm->NZ - 1;

			 static thread_local std::vector<real> a_diffySq(_iry * 2 + 1);
			 static thread_local std::vector<real> a_diffzSq(_irz * 2 + 1);
			 static thread_local std::vector<int> a_ipy(_iry * 2 + 1);
			 //int yc = 0; int zc = 0;
			 for (int l = ymin, yc = 0; l <= ymax; l++) {
				 const real my = l * pm->dy;
				 a_ipy[yc] = precalc_lat_y()[l + pm->NY];
				 const real diffy = p_diff_y(my, c->y());
				 a_diffySq[yc] = diffy*diffy;
				 yc++;
			 }

			 for (int m = zmin, zc = 0; m <= zmax; m++) {
				 const real mz = m * pm->dz;
				 const real diffz = mz - c->z();
				 a_diffzSq[zc] = diffz*diffz;
				 zc++;
			 }
			 //yc = 0; zc = 0;

			 for (int k = xmin; k <= xmax; k++) {
				 const real mx = k * pm->dx;
				 const int ipx = precalc_lat_x()[k + pm->NX];
				 const real diffx = p_diff_x(mx, c->x());
				 const real diffxSq = diffx*diffx;
				 for (int yc = 0; yc < _iry * 2 + 1; yc++) {
					 const int ipy = a_ipy[yc];
					 const real tmps = diffxSq + a_diffySq[yc];
					 for (int zc = 0, ipz = zmin; ipz <= zmax; ipz++, zc++) {
						 const real distSq = tmps + a_diffzSq[zc];
						 if (distSq < cradSq) {
							 //now for afm
							 cmap2.at(ipx,ipy,ipz) = 1;

							 //printf("%d %d %d\n", ipx, ipy, ipz);
							 if (distSq < normal_radSq) {
								 cmap1.at(ipx,ipy,ipz) = c;
							 }

						 }
					 }
				 }
			 }
		 }
	});
	tbb::parallel_for<size_t>(0, pm->NZ, [&cmap1, &cmap2](size_t l) {
        for (size_t j = 0; j<=pm->NX; j++) {
			cmap2.at(j,pm->NY,l) = cmap2.at(j,0,l);
			cmap1.at(j,pm->NY,l) = cmap1.at(j,0,l);
		}

        for (size_t k = 0; k <= pm->NY; k++) {
			cmap2.at(pm->NX,k,l) = cmap2.at(0,k,l);
			cmap1.at(pm->NX,k,l) = cmap1.at(0,k,l);
		}
	});
}
