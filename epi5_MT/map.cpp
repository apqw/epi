#include "map.h"
#include <tbb/parallel_for.h>
#include <tbb/blocked_range3d.h>
#include "cellmanager.h"
#include "cell.h"
#include "utils.h"


inline void set_lattice(Cell*const RESTRICT c) {
	c->lat[0] = precalc_lat_x()[(int)(c->x()*cont::inv_dx) + cont::NX];
	c->lat[1] = precalc_lat_y()[(int)(c->y()*cont::inv_dy) + cont::NY];
	c->lat[2] = precalc_lat_z()[(int)(c->z()*cont::inv_dz) + cont::NZ];
}

void setup_map_lat(CellManager & cman, FArr3D<const Cell*>& cmap1, FArr3D<cmask_ty>& cmap2)
{
	
	using namespace cont;
	tbb::parallel_for(tbb::blocked_range<int>(0, NX + 1), [&cmap1,&cmap2](const tbb::blocked_range< int >& range) {
		std::memset(&(cmap1()[range.begin()]), 0, sizeof(Cell*)*range.size()*(NY + 1)*(NZ + 1));
		std::memset(&(cmap2()[range.begin()]), 0, sizeof(cmask_ty)*range.size()*(NY + 1)*(NZ + 1));
	});
	cman.all_foreach_parallel_native([](Cell*const RESTRICT c) {
		set_lattice(c);
	});
	cman.all_foreach_parallel_native([&](const Cell*const c) {
		//auto&c = cman[i];
		
		auto& clat = c->lat;
		/*
		int  k, l, m, ipx, ipy, ipz;
		double mx, my, mz;
		*/
		const bool afm_cell = (c->state_mask()&(ALIVE_M | FIX_M | MUSUME_M))!=0u;
		
		 if (!afm_cell) {
			constexpr int _irx = irx;
			constexpr int _iry = iry;
			constexpr int _irz = irz;
			//double diffx, diffxSq, diffy, diffz, distSq;
			//const double crad = FAC_MAP*c->radius;
			//const double cradSq = crad*crad;
			const double normal_radSq = c->radius*c->radius;
			const int xmax = clat[0] + _irx; const int xmin = clat[0] - _irx;
			const int ymax = clat[1] + _iry; const int ymin = clat[1] - _iry;
			const int z_b = clat[2] - _irz;
			const int z_u = clat[2] + _irz;
			const int zmin = z_b >= 1 ? z_b : 1;
			const int zmax = z_u < (int)NZ ? z_u : (int)NZ - 1;

			static thread_local double a_diffySq_n[_iry * 2 + 1];
			static thread_local double a_diffzSq_n[_irz * 2 + 1];
			static thread_local int a_ipy_n[_iry * 2 + 1];
			//int yc = 0; int zc = 0;
			for (int l = ymin, yc = 0; l <= ymax; l++) {
				const double my = l * dy;
				a_ipy_n[yc] = precalc_lat_y()[l + NY];
				const double diffy = p_diff_y(my, c->y());
				a_diffySq_n[yc] = diffy*diffy;
				yc++;
			}

			for (int m = zmin, zc = 0; m <= zmax; m++) {
				const double mz = m * dz;
				const double diffz = mz - c->z();
				a_diffzSq_n[zc] = diffz*diffz;
				zc++;
			}
			//yc = 0; zc = 0;

			for (int k = xmin; k <= xmax; k++) {
				const double mx = k * dx;
				const int ipx = precalc_lat_x()[k + NX];
				const double diffx = p_diff_x(mx, c->x());
				const double diffxSq = diffx*diffx;
				for (int yc = 0; yc < _iry * 2 + 1; yc++) {
					const int ipy = a_ipy_n[yc];
					const double tmps = diffxSq + a_diffySq_n[yc];
					for (int zc = 0, ipz = zmin; ipz <= zmax; ipz++, zc++) {
						const double distSq = tmps + a_diffzSq_n[zc];
						if (distSq < normal_radSq) {
								cmap1()[ipx][ipy][ipz] = c;
						}
					}
				}
			}
		}else {
			 constexpr int _irx = 2 * irx;
			 constexpr int _iry = 2 * iry;
			 constexpr int _irz = 2 * irz;
			 //double diffx, diffxSq, diffy, diffz, distSq;
			 const double crad = FAC_MAP*c->radius;
			 const double cradSq = crad*crad;
			 const double normal_radSq = c->radius*c->radius;
			 const int xmax = clat[0] + _irx; const int xmin = clat[0] - _irx;
			 const int ymax = clat[1] + _iry; const int ymin = clat[1] - _iry;
			 const int z_b = clat[2] - _irz;
			 const int z_u = clat[2] + _irz;
			 const int zmin = z_b >= 1 ? z_b : 1;
			 const int zmax = z_u < (int)NZ ? z_u : (int)NZ - 1;

			 static thread_local double a_diffySq[_iry * 2 + 1];
			 static thread_local double a_diffzSq[_irz * 2 + 1];
			 static thread_local int a_ipy[_iry * 2 + 1];
			 //int yc = 0; int zc = 0;
			 for (int l = ymin, yc = 0; l <= ymax; l++) {
				 const double my = l * dy;
				 a_ipy[yc] = precalc_lat_y()[l + NY];
				 const double diffy = p_diff_y(my, c->y());
				 a_diffySq[yc] = diffy*diffy;
				 yc++;
			 }

			 for (int m = zmin, zc = 0; m <= zmax; m++) {
				 const double mz = m * dz;
				 const double diffz = mz - c->z();
				 a_diffzSq[zc] = diffz*diffz;
				 zc++;
			 }
			 //yc = 0; zc = 0;

			 for (int k = xmin; k <= xmax; k++) {
				 const double mx = k * dx;
				 const int ipx = precalc_lat_x()[k + NX];
				 const double diffx = p_diff_x(mx, c->x());
				 const double diffxSq = diffx*diffx;
				 for (int yc = 0; yc < _iry * 2 + 1; yc++) {
					 const int ipy = a_ipy[yc];
					 const double tmps = diffxSq + a_diffySq[yc];
					 for (int zc = 0, ipz = zmin; ipz <= zmax; ipz++, zc++) {
						 const double distSq = tmps + a_diffzSq[zc];
						 if (distSq < cradSq) {
							 //now for afm
							 cmap2()[ipx][ipy][ipz] = 1;

							 //printf("%d %d %d\n", ipx, ipy, ipz);
							 if (distSq < normal_radSq) {
								 cmap1()[ipx][ipy][ipz] = c;
							 }

						 }
					 }
				 }
			 }
		 }
	});
	tbb::parallel_for<size_t>(0, NZ, [&cmap1, &cmap2](size_t l) {
        for (size_t j = 0; j<=NX; j++) {
			cmap2()[j][NY][l] = cmap2()[j][0][l];
			cmap1()[j][NY][l] = cmap1()[j][0][l];
		}

        for (size_t k = 0; k <= NY; k++) {
			cmap2()[NX][k][l] = cmap2()[0][k][l];
			cmap1()[NX][k][l] = cmap1()[0][k][l];
		}
	});
}
