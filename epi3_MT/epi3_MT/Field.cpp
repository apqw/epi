#define _CRT_SECURE_NO_WARNINGS

#include "Field.h"
#include <string>
#include <cassert>
#include <iostream>
void Field::interact_cell() {

	
	cells.foreach_parallel_native([](CellPtr& c) {
			switch (c->state()) {
			case MEMB:
				if (c->pos[2]() < c->radius()) {
					c->wall_interact();
				}
				c->MEMB_interact();
				break;
			case DER:
				if (c->pos[2]() < c->radius()) {
					c->wall_interact();
				}
				c->DER_interact();
				break;
			case ALIVE:case AIR:case DEAD:
				c->AL_AIR_DE_interact();
				break;
			case FIX:
				c->FIX_interact();
				break;
			case MUSUME:
				c->MUSUME_interact();
				
				break;
			default:
				break;
			}
			c->pair_interact();
		});
	
}
//do not multithread
void Field::cell_state_renew() {
	cells.other_foreach([&](CellPtr& c, int i) {
		switch (c->state())
		{
        case MUSUME:
            c->MUSUME_state_renew();
            break;
		case FIX:
			c->FIX_state_renew();
			break;

		case DEAD:case AIR:
			c->DEAD_AIR_state_renew();
			break;
		case ALIVE:
			c->ALIVE_state_renew();
		default:
			break;
		}
		if (c->should_deleted()) {
			printf("remove detected\n");
			cells.remove_queue(i);
		}
		if (c->has_new_pair()) {
			printf("new pair detected\n");
			cells.add_queue(c->pair);
			c->pair->pair = c;
			c->set_as_no_more_new_pair();
		}
	
	});
}
/*
	fix next value
*/
void Field::cell_pos_periodic_fix() {
	using namespace cont;
	cells.foreach_parallel_native([](CellPtr& c) {
		if (c->pos[0].get_next_value() > LX) {
			c->pos[0].force_set_next_value(c->pos[0].get_next_value() - LX);
		}
		else if (c->pos[0].get_next_value() < 0) {
			c->pos[0].force_set_next_value(c->pos[0].get_next_value() + LX);
		}

		if (c->pos[1].get_next_value() > LY) {
			c->pos[1].force_set_next_value(c->pos[1].get_next_value() - LY);
		}
		else if (c->pos[1].get_next_value() < 0) {
			c->pos[1].force_set_next_value(c->pos[1].get_next_value() + LY);
		}
	});
}

void Field::connect_cells() {
	using namespace cont;
    static RawArr3D<std::atomic<int>,ANX,ANY,ANZ> aindx = { };
	static RawArr3D<RawArr1D<Cell*, N3>, ANX, ANY, ANZ> area = { nullptr };
	static bool mflg = false;
	for (int i = 0; i < ANX; i++) {
		for (int j = 0; j < ANY; j++) {
			for (int k = 0; k < ANZ; k++) {
				aindx[i][j][k] = 0;
			}
		}
	}
	cells.foreach_parallel_native([&](CellPtr& c) {
		int aix, aiy, aiz;
		aix = (int)((0.5*LX - p_diff_sc_x(0.5*LX, c->pos[0]())) / AREA_GRID);
		aiy = (int)((0.5*LY - p_diff_sc_y(0.5*LY, c->pos[1]())) / AREA_GRID);
		aiz = (int)(min0(c->pos[2]()) / AREA_GRID);

		//assert();
		if ((aix >= ANX || aiy >= ANY || aiz >= ANZ || aix < 0 || aiy < 0 || aiz < 0)) {
			printf("err\n");
			printf("cx:%lf cy:%lf cz:%lf\n",c->pos[0](), c->pos[1](), c->pos[2]());
			printf("aix:%d aiy:%d aiz:%d\n", aix, aiy, aiz);
			assert(false);
		}
		int atm = aindx[aix][aiy][aiz]++;
		area[aix][aiy][aiz][atm] = c.get();
		
		assert(aindx[aix][aiy][aiz] < N3);

	});
	cells.foreach_parallel_native([&](CellPtr& c) {
			c->connected_cell.set_count(c->state()==MEMB?4:0);
			int anx = (int)(c->pos[0]() / AREA_GRID);
			int any = (int)(c->pos[1]() / AREA_GRID);
			int anz = (int)(c->pos[2]() / AREA_GRID);
			
			assert(!(anx >= ANX || any >= ANY || anz >= ANZ || anx < 0 || any < 0 || anz < 0));
			Vec3<double> diffv;
			double rad_sum; double diffx, diffy, diffz;
			int ii = 2, aix, aiy, aiz;
			for (int j = anx - ii; j <= anx + ii; j++) {
				aix = j;
				if (j < 0) {
					aix += ANX;
				}
				else if (aix >= ANX) {
					aix -= ANX;
				}

				for (int k = any - ii; k <= any + ii; k++) {
					aiy = k;
					if (k < 0) {
						aiy += ANY;
					}
					else if (k >= ANY) {
						aiy -= ANY;
					}

					for (int l = anz - ii; l <= anz + ii; l++) {
						aiz= l;
						if (l < 0) {
							aiz += ANZ;
						}
						else if (l >= ANZ) {
							aiz -= ANZ;
						}

						for (int m = 0; m < aindx[aix][aiy][aiz]; m++) {
							Cell* o = area[aix][aiy][aiz][m];
							if (c.get() == o||(c->state()==MEMB && o->state()==MEMB))continue;
						
							//diffv = c->pos - o->pos;
							diffx = p_diff_sc_x(c->pos[0]() , o->pos[0]());
							diffy = p_diff_sc_y(c->pos[1](), o->pos[1]());
							diffz=c->pos[2]()-o->pos[2]();
							rad_sum = c->radius() + o->radius();
							/*
							if (fabs(tmp - tmp2) > 1) {
								printf("%lf\n", tmp - tmp2);
							}
							*/
							if (diffx*diffx + diffy*diffy + diffz*diffz <= LJ_THRESH*LJ_THRESH*rad_sum*rad_sum) {
								//printf("connecting... %d \n", c->connected_cell.count()+1);
								c->connected_cell.add(o);
								assert(c->connected_cell.count() < N2);
							}
							
						}
					}
				}
			}

	});
	cells.foreach_parallel_native([&](CellPtr& c) {

		//delete unconnected cell value
		for (auto it = c->gj.begin(); it != c->gj.end();) {
			if (!c->connected_cell.exist(it->first)) {
				it = c->gj.erase(it);
			}
			else {
				++it;
			}
		}

		//set newly connected cell value
		c->connected_cell.foreach([&](Cell* cptr) {
			if (c->gj.count(cptr)==0) {
				c->gj.emplace(cptr, w0);
			}
		});
	
	});
}

void Field::set_cell_lattice()
{
	cells.foreach_parallel_native([](CellPtr& c) {
		c->set_lattice();
	});
}

double Field::calc_zzmax()
{
	double zmax = 0;
	cells.other_foreach([&zmax](CellPtr& c, int i) {
		auto& st = c->state();
		if (st == DEAD || st == ALIVE || st == MUSUME || st == FIX) {
			if (zmax < c->pos[2]())zmax = c->pos[2]();
		}
	});
	return zmax;
}

void Field::cell_dynamics() {
	cells.memb_foreach_parallel_native([](CellPtr& memb) {
		memb->memb_bend_calc1();
	});

	cells.memb_foreach_parallel_native([](CellPtr& memb) {
		memb->memb_bend_calc2();
	});

	cells.memb_foreach_parallel_native([](CellPtr& memb) {
		memb->memb_bend_interact();
	});

	interact_cell();
	cells.all_cell_update();

	cell_state_renew();
	cells.update();
	cells.all_cell_update();

	cells.other_foreach([](CellPtr& c, int i) {
        if (c->pair != nullptr&&c>c->pair) {
			c->pair_disperse();

		}
	});

	cell_pos_periodic_fix();
	cells.all_cell_update();
	connect_cells();
}

void Field::main_loop()
{
	for (int i = 0; i < cont::NUM_ITR; i++) {
        if(i%100==0)printf("loop:%d\n", i);
		cells.all_cell_update();
		cell_dynamics();
		zzmax = calc_zzmax();
		set_cell_lattice();
		setup_map();
		calc_b();
		b_update();
		if (i % 10000 == 0) {
			printf("calc ca...\n");
			calc_ca();
			
			printf("end.\n");
		}
	}
}

void Field::setup_map()
{
	using namespace cont;

		for (int i = 0; i != NX; ++i) {
			for (int j = 0; j <= NY; j++) {
				for (int k = 0; k <= NZ; k++) {
					air_stim_flg[i][j][k] = 0;
					cell_map[i][j][k] = nullptr;
					cell_map2[i][j][k] = 0;
				}
			}
		}
	
	cells.foreach_parallel_native([&](CellPtr& c) {
		
		auto& cv = c->pos;
		//double crad = c->old_data.radius;
		auto& clat = c->lat;
		int  k, l, m, ipx, ipy, ipz;
		double mx, my, mz;

		constexpr int _irx = 2 * irx;
		constexpr int _iry = 2 * iry;
		constexpr int _irz = 2 * irz;
		double diffx,diffxSq, diffy, diffz, distSq;
		double crad = FAC_MAP*c->radius();
		double cradSq = crad*crad;
		double normal_radSq = c->radius()*c->radius();
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
			ipy = l;
			if (l < 0) {
				ipy += NY;
			}
			else if (l >= NY) {
				ipy -= NY;
			}
			a_ipy[yc] = ipy;
			diffy = my - cv[1]();
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
			diffz = mz - cv[2]();
			a_diffzSq[zc] = diffz*diffz;
			zc++;
		}
		yc = 0; zc = 0;
		for (k = xmin; k <= xmax; k++) {
			//imx = k;
			mx = k * dx;

			ipx = k;


			if (k < 0) {
				ipx += NX;
			}
			else if (k >= NY) {
				ipx -= NX;
			}

			diffx = mx - cv[0]();
			if (fabs(diffx) >= 0.5*LX) {
				if (diffx > 0) {
					diffx -= LX;
				}
				else {
					diffx += LX;
				}
			}
			diffxSq = diffx*diffx;
			for (yc=0,l = ymin; l <= ymax; l++,yc++) {
				ipy = a_ipy[yc];

				
				for (zc=0,m = zmin; m <= zmax; m++,zc++) {
					//ipz = imz = iz + m;
					ipz = m;
					if ((distSq = diffxSq + a_diffySq[yc] + a_diffzSq[zc]) < cradSq) {
						cell_map2[ipx][ipy][ipz] = 1;

						if (distSq < normal_radSq) {
							cell_map[ipx][ipy][ipz] = c.get();
							air_stim_flg[ipx][ipy][ipz] = (c->state() == AIR);
							
							
						}
					}
				}
			}
		}

	});
}

void Field::calc_b() {
	using namespace cont;
	int iz_bound = (int)((zzmax + FAC_MAP*R_max) / dz);

	tbb::parallel_for(tbb::blocked_range<int>(0, NX ), [&](const tbb::blocked_range< int >& range) {
		for (int j = range.begin(); j != range.end(); ++j) {
			int prev_x = j - 1;
			int next_x = j + 1;
			if (prev_x < 0) {
				prev_x += NX;
			}
			if (next_x >= NX) {
				next_x -= NX;
			}
			for (int k = 0; k < NY; k++) {
				int prev_y = k - 1;
				int next_y = k + 1;
				if (prev_y < 0) {
					prev_y += NY;
				}
				if (next_y >= NY) {
					next_y -= NY;
				}
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
					double dum_age = 0;
					bool flg_cornified = false;
					if (cell_map[j][k][l]!=nullptr) {
						if (get_state_mask(cell_map[j][k][l] ->state())&(ALIVE_M|DEAD_M)) {
							dum_age = cell_map[j][k][l]->agek();
							flg_cornified = true;
						}
					}
					ext_stim[j][k][l]+=DT_Ca*(DB * (cell_map2[prev_x][k][l] * (ext_stim[prev_x][k][l]() - ext_stim[j][k][l]())
						+ cell_map2[j][prev_y][l] * (ext_stim[j][prev_y][l]() - ext_stim[j][k][l]())
						+ cell_map2[j][k][prev_z] * (ext_stim[j][k][prev_z]() - ext_stim[j][k][l]())
						+ cell_map2[next_x][k][l] * (-ext_stim[j][k][l]() + ext_stim[next_x][k][l]())
						+ cell_map2[j][next_y][l] * (-ext_stim[j][k][l]() + ext_stim[j][next_y][l]())
						+ cell_map2[j][k][next_z] * (-ext_stim[j][k][l]() + ext_stim[j][k][next_z]())) / (dz * dz)
						+ fB(dum_age, ext_stim[j][k][l]() ,flg_cornified));


				}
			}
		}
	});
	for (int l = 0; l <= iz_bound; l++) {
		for (int j = 0; j < NX; j++) ext_stim[j][NY][l].force_set_next_value(ext_stim[j][0][l].get_next_value());
		for (int k = 0; k <= NY; k++) ext_stim[NX][k][l].force_set_next_value(ext_stim[0][k][l].get_next_value());
	}
}

void Field::b_update()
{
	for (auto& x : ext_stim) {
		for (auto& y : x) {
			for (auto& z : y) {
				z.update();
			}
		}
	}
}
double th(CELL_STATE state, double age) {
	assert(get_state_mask(state)&(ALIVE_M | FIX_M | MUSUME_M));
	using namespace cont;
	if (state == ALIVE) {
		return thgra + ((thpri - thgra)*0.5) * (1.0 + tanh((THRESH_SP - age) / delta_th));
	}
	else if(get_state_mask(state)&( FIX_M | MUSUME_M)) {
		return thpri;
	}
	return 0;
}
void Field::calc_ca()
{
	cells.other_foreach([](CellPtr& c, int i) {
		if (get_state_mask(c->state())&(ALIVE_M|FIX_M|MUSUME_M)) {
			c->ca2p_avg.force_set_next_value(0);
		}
	});
	using namespace cont;
	int iz_bound = (int)((zzmax + FAC_MAP*R_max) / dz);
	for (int cstp = 0; cstp <Ca_ITR; cstp++) {

		cells.other_foreach_parallel_native([&iz_bound,this](CellPtr& c) {
			if (c->state() == DEAD) {
				
				int count = 0;
				double tmp = 0;
				c->connected_cell.foreach([&count,&c,&tmp](Cell* conn) {
					if (conn->state() == ALIVE) {
						tmp += conn->IP3();
						count++;
					}
				});
				tmp = DT_Ca*(dp*(tmp-count*c->IP3())-Kpp*c->IP3());
				c->IP3 += tmp;
			}
			else if (get_state_mask(c->state())&(ALIVE_M | FIX_M | MUSUME_M)) {
				int ix = c->lat[0]; int iy = c->lat[1]; int iz = c->lat[2];

				assert(iz < iz_bound);

				double tmp_a = grid_avg8(ATP, ix, iy, iz);
				double tmp_B = grid_avg8(ext_stim, ix, iy, iz);

				c->diffu = fu(c->ca2p(), c->ex_inert(), c->IP3(), tmp_B);
				double _th = thpri;
				double _Kpa = Kpri;
				double IAGv = iage_kitei;
				if (c->state()==ALIVE){
					//ALIVE
					_th = thgra + ((thpri - thgra)*0.5) * (1.0 + tanh((THRESH_SP - c->agek()) / delta_th));
					_Kpa = Kgra + ((Kpri - Kgra)*0.5) * (1.0 + tanh((THRESH_SP - c->agek()) / delta_K));
					IAGv= 0.5*(1.0 + tanh((c->agek() - THRESH_SP) / delta_I));
				}
				c->ex_inert += DT_Ca*((para_k2*para_k2 / (para_k2*para_k2 + c->ca2p()*c->ca2p()) - c->ex_inert()) / _th);
				c->IP3 += DT_Ca*(_Kpa*tmp_a / (H0 + tmp_a) - Kpp*c->IP3());
				double tmp_diffu = 0;
				double tmp_IP3 = 0;
				c->connected_cell.foreach([&c,&tmp_diffu,&tmp_IP3,IAGv](Cell* conn) {
					if (get_state_mask(conn->state())&(ALIVE_M | DEAD_M | FIX_M | MUSUME_M)) {
						tmp_diffu += c->gj.at(conn)()*(conn->ca2p() - c->ca2p());
						tmp_IP3 += c->gj.at(conn)()*(conn->IP3() - c->IP3());
						if (conn->state() == ALIVE) {
							c->gj.at(conn) += DT_Ca*fw(fabs(conn->ca2p() - c->ca2p()), c->gj.at(conn)());
						}
					}
				});
				c->diffu+=ca2p_du*IAGv*tmp_diffu;
				c->IP3+= DT_Ca*dp*IAGv*tmp_IP3;
				c->ca2p += DT_Ca*c->diffu;
				c->ca2p_avg += c->ca2p()+ DT_Ca*c->diffu;
			}
		});

		tbb::parallel_for(tbb::blocked_range<int>(0, NX), [&](const tbb::blocked_range< int >& range) {
			for (int j = range.begin(); j!= range.end(); ++j) {
				int prev_x = j - 1;
				if (prev_x < 0) {
					prev_x += NX;
				}
				int next_x = j + 1;
				if (next_x >= NX) {
					next_x -= NX;
				}
				for (int k = 0; k < NY; k++) {
					int prev_y = k - 1;
					if (prev_y < 0) {
						prev_y += NY;
					}
					int next_y = k + 1;
					if (next_y >= NY) {
						next_y -= NY;
					}
					for (int l = 0; l < iz_bound; l++) {
						int prev_z = l - 1;
						if (l == 0)prev_z = 1;
						int next_z = l + 1;
						if (l == NZ)next_z = NZ - 1;
						double 	tmp = 0;
						if (cell_map[j][k][l] != nullptr) {
							if (get_state_mask(cell_map[j][k][l]->state())&(ALIVE_M | FIX_M | MUSUME_M)) {
								tmp = cell_map[j][k][l]->diffu;
							}
						}
						cell_map[j][k][l] != nullptr ? cell_map[j][k][l]->diffu : 0.0;
						ATP[j][k][l]
							+= DT_Ca*(Da * (cell_map2[prev_x][k][l] * (ATP[prev_x][k][l] - ATP[j][k][l])
								+ cell_map2[next_x][k][l] * (ATP[next_x][k][l] - ATP[j][k][l])
								+ cell_map2[j][prev_y][l] * (ATP[j][prev_y][l] - ATP[j][k][l])
								+ cell_map2[j][next_y][l] * (ATP[j][next_y][l] - ATP[j][k][l])
								+ cell_map2[j][k][prev_z] * (ATP[j][k][prev_z] - ATP[j][k][l])
								+ cell_map2[j][k][next_z] * (ATP[j][k][next_z] - ATP[j][k][l])) / (dx*dx)
								+ fa(tmp, ATP[j][k][l]()));
						ATP[j][k][l] += DT_Ca*air_stim_flg[j][k][l] * AIR_STIM;
					}
				}
			}
		});
		for (int l = 0; l <= iz_bound; l++) {
			for (int j = 0; j < NX; j++) ATP[j][NY][l].force_set_next_value(ATP[j][0][l].get_next_value());
			for (int k = 0; k <= NY; k++) ATP[NX][k][l].force_set_next_value(ATP[0][k][l].get_next_value());
		}
		cells.other_foreach_parallel_native([&iz_bound](CellPtr& c) {
			c->ca2p.update();
			c->ex_inert.update();
			c->IP3.update();
			c->diffu = 0;
			for (auto& gjv : c->gj) {
				gjv.second.update();
			}
		});
		ATP_update();
	}
	cells.other_foreach_parallel_native([](CellPtr& c) {
		if (get_state_mask(c->state())&(ALIVE_M | FIX_M | MUSUME_M)) {
			c->ca2p_avg /= Ca_ITR;
		}
		else {
			c->ca2p_avg.force_set_next_value(0);
			c->ca2p.force_set_next_value(0);
		}
	});
}

void Field::ATP_update()
{
	for (auto& x : ATP) {
		for (auto& y : x) {
			for (auto& z : y) {
				z.update();
			}
		}
	}
}

void Field::init_with_file(std::ifstream& dstrm) {
	using namespace cont;
	std::string line;
	CELL_STATE state;
	int div_times, touch, pair_cell_id, stem_orig_id;
	double rad, ageb, agek, x, y, z, fat, spr_len, ex_fat;
	int id_count = 0;
	std::vector<CELL_STATE> tmp_pair_state(MAX_CELL_NUM, UNUSED);
	std::vector<CellPtr> tmp_pair(MAX_CELL_NUM, nullptr);
	unsigned int phase = 0;
	int nmemb = 0;
	int nder = 0;
	while (!dstrm.eof()) {
		std::getline(dstrm, line);
        sscanf(line.c_str(), "%*d %d %lf %lf %lf %*f %lf %lf %lf %*f %d %lf %lf %d %lf %d %d",
			&state, &rad, &ageb, &agek, &x, &y, &z, &div_times, &ex_fat, &fat, &touch, &spr_len, &pair_cell_id, &stem_orig_id);
		if (state == BLANK)break; //owari
		if (SYSTEM == BASAL && (state == ALIVE || state == DEAD || state == AIR)) {
			printf(" input date must not contain ALIVE or DEAD in case of BASAL\n");
			exit(1);
		}
		if (state == DER && rad != R_der) {
			printf("radii of DER not consistent with param.h\n");
			exit(1);
		}
		if (state == MEMB && rad != R_memb) {
			printf("radii of DER not consistent with param.h\n");
			exit(1);
		}
		if (phase == 0 && state != MEMB) {
			assert(state == DER);
			phase++;
		}

		if (phase == 1 && state != DER) {
			//assert(state == DER);
			phase++;
		}
		if (phase > 0 && state == MEMB) {
			printf("non phase0 memb\n");
			exit(1);
		}
		if (phase > 1 && state == DER) {
			printf("non phase1 der\n");
			exit(1);
		}

		//if (state == FIX)printf("FIX\n");
		printf("Phase %d  Cell loaded:%d\n",phase, id_count++);
		if (state == MEMB)nmemb++;
		if (state == DER)nder++;
		auto cptr = std::make_shared<Cell>(
			state,
			std::initializer_list<DV<double>>{x, y, z},
			rad, //radius
			u0, //ca2p 
			u0, //ca2pavg (initial value unused)
			p0, //IP3
			v0, //ex_inert
            agek,
            #ifdef AGE_DBG
                    5,
            #else
                    ageb,
            #endif
			ex_fat, fat,//ex_fat in_fat
			spr_len,
			0,//spring force (unused)
            state == FIX ? agki_max_fix :
			state == MUSUME ?
#ifdef AGE_DBG
			agki_max*0.01 :
#else
		agki_max:
#endif
			0,//div thresh
			0,//poisson??
			state==FIX?div_max:0,//rest div times
			stem_orig_id<MALIG_NUM,//malignant
			touch == 1//touch
			);
		cells.add_direct(cptr);

		if (pair_cell_id > -1) {
			if (tmp_pair[id_count] != nullptr) {
				assert(tmp_pair[id_count]->pair == nullptr);
				cptr->pair = tmp_pair[id_count];
				tmp_pair[id_count]->pair = cptr;
				/*
				if (state == FIX || tmp_pair_state[id_count] == FIX) {
					cptr->old_data.spring_force =
						tmp_pair[id_count]->old_data.spring_force = Kspring;
				}
				else if (state == MUSUME &&  tmp_pair_state[id_count] == MUSUME) {
					cptr->old_data.spring_force =
						tmp_pair[id_count]->old_data.spring_force = Kspring_d;
				}
				else {
					cptr->old_data.spring_force =
						tmp_pair[id_count]->old_data.spring_force = 0;
				}
				*/
			}
			else {
				tmp_pair[pair_cell_id] = cptr;
				tmp_pair_state[pair_cell_id] = state;
			}
		}


	}

	cells.set_der_num(nder);
	cells.set_memb_num(nmemb);
	using namespace cont;
	//memb indices init
	int  jj, kk;
	auto& raw_cells = cells._raw_cell_set();
	printf("setting memb initial connection...\n");
	for (int j = 0; j < nmemb; j++) {
		auto& cptr = raw_cells[j];
		jj = j%NMX;
		kk = j / NMX;
		if (jj == 0) {
			cptr->mbd.memb_l = raw_cells[j + NMX - 1];
		}
		else {
			cptr->mbd.memb_l = raw_cells[j - 1];
		}

		if (jj <= 1) {
			cptr->mbd.memb_ll = raw_cells[j + NMX - 2];
		}
		else {
			cptr->mbd.memb_ll = raw_cells[j - 2];
		}

		if (jj == NMX-1) {
			cptr->mbd.memb_r = raw_cells[j - (NMX - 1)];
		}
		else {
			cptr->mbd.memb_r = raw_cells[j + 1];
		}

		if (kk == 0) {
			cptr->mbd.memb_b= raw_cells[j +NMX*NMY-NMX];
		}
		else {
			cptr->mbd.memb_b = raw_cells[j -NMX];
		}

		if (kk <= 1) {
			cptr->mbd.memb_bb = raw_cells[j + NMX*NMY - 2*NMX];
		}
		else {
			cptr->mbd.memb_bb = raw_cells[j - 2*NMX];
		}

		if (kk == NMY-1) {
			cptr->mbd.memb_u = raw_cells[j -(NMX*NMY -  NMX)];
		}
		else {
			cptr->mbd.memb_u = raw_cells[j + NMX];
		}

		cptr->connected_cell.add(cptr->mbd.memb_l);
		cptr->connected_cell.add(cptr->mbd.memb_r);
		cptr->connected_cell.add(cptr->mbd.memb_b);
		cptr->connected_cell.add(cptr->mbd.memb_u);
	}
	printf("done.\n");
	printf("setting initial connection...\n");
	connect_cells();
	printf("done.\n");
}
