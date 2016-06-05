#include "map.h"
#include <tbb/parallel_for.h>
#include <tbb/blocked_range3d.h>
#include "cellmanager.h"
#include "cell.h"
#include "utils.h"



void set_lattice(Cell* c) {
	c->lat[0] = precalc_lat_x()[(int)(c->x()*cont::inv_dx) + cont::NX];
	c->lat[1] = precalc_lat_y()[(int)(c->y()*cont::inv_dy) + cont::NY];
	c->lat[2] = precalc_lat_z()[(int)(c->z()*cont::inv_dz) + cont::NZ];
}

void setup_map_lat(CellManager & cman, Field<Cell*,cont::NX + 1, cont::NY + 1, cont::NZ + 1>& cmap1, Field<uint_fast8_t,cont::NX + 1, cont::NY + 1, cont::NZ + 1>& cmap2)
{

	using namespace cont;
	tbb::parallel_for(tbb::blocked_range<int>(0, NX + 1), [&](const tbb::blocked_range< int >& range) {
		std::memset(&(cmap1()[range.begin()]), 0, sizeof(Cell*)*range.size()*(NY + 1)*(NZ + 1));
		std::memset(&(cmap2()[range.begin()]), 0, sizeof(uint_fast8_t)*range.size()*(NY + 1)*(NZ + 1));
	});

	cman.all_foreach_parallel([&](size_t i) {
		auto&c = cman[i];
		set_lattice(c);
		auto& clat = c->lat;
		int  k, l, m, ipx, ipy, ipz;
		double mx, my, mz;

		constexpr int _irx = 2 * irx;
		constexpr int _iry = 2 * iry;
		constexpr int _irz = 2 * irz;
		double diffx, diffxSq, diffy, diffz, distSq;
		double crad = FAC_MAP*c->radius;
		double cradSq = crad*crad;
		double normal_radSq = c->radius*c->radius;
		int xmax = clat[0] + _irx; int xmin = clat[0] - _irx;
		int ymax = clat[1] + _iry; int ymin = clat[1] - _iry;
		int z_b = clat[2] - _irz;
		int z_u = clat[2] + _irz;
		int zmin = z_b >= 1 ? z_b : 1;
		int zmax = z_u < NZ ? z_u : NZ - 1;

		double a_diffySq[_iry * 2 + 1];
		double a_diffzSq[_irz * 2 + 1];
		int a_ipy[_iry * 2 + 1];
		int yc = 0; int zc = 0;
		for (l = ymin; l <= ymax; l++) {
			my = l * dy;
			a_ipy[yc] = precalc_lat_y()[l+NY];
			diffy = my - c->y();
			if (fabs(diffy) >= 0.5*LY) {
				if (diffy > 0) {
					diffy -= LY;
				}
				else {
					diffy += LY;
				}
			}
			a_diffySq[yc] = diffy*diffy;
			yc++;
		}

		for (m = zmin; m <= zmax; m++) {
			mz = m * dz;
			diffz = mz - c->z();
			a_diffzSq[zc] = diffz*diffz;
			zc++;
		}
		yc = 0; zc = 0;
		bool afm_cell = c->state == ALIVE || c->state == FIX || c->state == MUSUME;
		for (k = xmin; k <= xmax; k++) {
			mx = k * dx;
			ipx = precalc_lat_x()[k+NX];
			diffx = mx - c->x();
			if (fabs(diffx) >= 0.5*LX) {
				if (diffx > 0) {
					diffx -= LX;
				}
				else {
					diffx += LX;
				}
			}
			diffxSq = diffx*diffx;
			for (yc = 0, l = ymin; l <= ymax; l++, yc++) {
				ipy = a_ipy[yc];
				for (zc = 0, ipz = zmin; ipz <= zmax; ipz++, zc++) {
					
					if ((distSq = diffxSq + a_diffySq[yc] + a_diffzSq[zc]) < cradSq) {
						if (afm_cell) {
							cmap2()[ipx][ipy][ipz] = 1;
						}
						
						if (distSq < normal_radSq) {
							cmap1()[ipx][ipy][ipz] = c;
						}
						
					}
				}
			}
		}
	});
	for (int l = 0; l<NZ; l++) {
		for (int j = 0; j<NX; j++) {
			cmap2()[j][NY][l] = cmap2()[j][0][l];
			cmap1()[j][NY][l] = cmap1()[j][0][l];
		}

		for (int k = 0; k <= NY; k++) {
			cmap2()[NX][k][l] = cmap2()[0][k][l];
			cmap1()[NX][k][l] = cmap1()[0][k][l];
		}
	}
}
