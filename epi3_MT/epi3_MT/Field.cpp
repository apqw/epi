#include "Field.h"
#include <string>
#include <cassert>
#include <iostream>
void Field::interact_cell() {

		cells.foreach_parallel([](CellPtr& c, int i) {
			switch (c->state()) {
			case MEMB:
				
				c->MEMB_interact();
				break;
			case DER:
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
		case FIX:
			c->FIX_state_renew();
			break;
		case MUSUME:
			c->MUSUME_state_renew();
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
			cells.remove_queue(i);
		}
		if (c->has_new_pair()) {
			printf("new pair detected\n");
			cells.add_queue(c->pair);
			c->pair->pair = c;
		}
	
	});
}
/*
	fix next value
*/
void Field::cell_pos_periodic_fix() {
	using namespace cont;
	cells.foreach_parallel([](CellPtr& c, int i) {
		if (c->pos[0] > LX) {
			c->pos[0].force_set_next_value(c->pos[0]() - LX);
		}
		else if (c->pos[0] < 0) {
			c->pos[0].force_set_next_value(c->pos[0]() + LX);
		}

		if (c->pos[1] > LY) {
			c->pos[1].force_set_next_value(c->pos[1]() - LY);
		}
		else if (c->pos[1] < 0) {
			c->pos[1].force_set_next_value(c->pos[1]() + LY);
		}
	});
}

void Field::connect_cells() {
	using namespace cont;
	static RawArr3D<std::atomic<int>,ANX,ANY,ANZ> aindx = { 0 };
	static RawArr3D<RawArr1D<Cell*, N3>, ANX, ANY, ANZ> area = { nullptr };
	static bool mflg = false;

	for (int i = 0; i < ANX; i++) {
		for (int j = 0; j < ANY; j++) {
			for (int k = 0; k < ANZ; k++) {
				aindx[i][j][k] = 0;
			}
		}
	}
	cells.foreach_parallel([&](CellPtr& c, int i) {
		int aix, aiy, aiz;
		aix = (int)((0.5*LX - p_diff_sc_x(0.5*LX, c->pos[0]())) / AREA_GRID);
		aiy = (int)((0.5*LY - p_diff_sc_y(0.5*LY, c->pos[1]())) / AREA_GRID);
		aiz = (int)(min0(c->pos[2]()) / AREA_GRID);

		assert(!(aix >= ANX || aiy >= ANY || aiz >= ANZ || aix < 0 || aiy < 0 || aiz < 0));
		int atm = aindx[aix][aiy][aiz]++;
		area[aix][aiy][aiz][atm] = c.get();
		
		assert(aindx[aix][aiy][aiz] < N3);

	});

	cells.foreach_parallel([&](CellPtr& c, int i) {
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
							if (c.get() <= o||(c->state()==MEMB && o->state()==MEMB))continue;
							//avoid double counting?
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
								o->connected_cell.add(c.get());
								assert(c->connected_cell.count() < N2);
							}
							
						}
					}
				}
			}

	});

}

void Field::cell_dynamics() {
	cells.memb_foreach_parallel([](CellPtr& memb, int i) {
		memb->memb_bend_calc1();
	});

	cells.memb_foreach_parallel([](CellPtr& memb, int i) {
		memb->memb_bend_calc2();
	});

	cells.memb_foreach_parallel([](CellPtr& memb, int i) {
		memb->memb_bend_interact();
	});

	interact_cell();
	cells.all_cell_update();

	cell_state_renew();
	cells.update();
	cells.all_cell_update();

	cells.other_foreach([](CellPtr& c, int i) {
		if (c->pair != nullptr) {
			c->pair_disperse();
			c->update();
		}
	});

	cell_pos_periodic_fix();
	cells.all_cell_update();
	connect_cells();
}

void Field::main_loop()
{
	for (int i = 0; i < cont::NUM_ITR; i++) {
		printf("loop:%d\n", i);
		cell_dynamics();
		//add lattice calc
		setup_map();
	}
}

void Field::setup_map()
{
	using namespace cont;
	printf("");
	
	tbb::parallel_for(tbb::blocked_range<int>(0,NX+1), [&](const tbb::blocked_range< int >& range) {
		for (int i = range.begin(); i != range.end(); ++i) {
			for (int j = 0; j <= NY; j++) {
				for (int k = 0; k <= NZ; k++) {
					age_map[i][j][k] = 0;
					cornif_map[i][j][k] = 0;
					air_stim_flg[i][j][k] = 0;
					cell_map[i][j][k] = nullptr;
					cell_map2[i][j][k] = 0;
				}
			}
		}
	});
	cells.foreach_parallel([&](CellPtr& c, int i) {
		
		auto& cv = c->pos;
		//double crad = c->old_data.radius;
		auto& clat = c->lat;
		int  k, l, m, ipx, ipy, ipz, imx, imy;
		double mx, my, mz;

		int _irx = 2 * irx;
		int _iry = 2 * iry;
		int _irz = 2 * irz;
		double diffx, diffy, diffz, distSq;
		double crad = FAC_MAP*c->radius();
		double cradSq = crad*crad;
		double normal_radSq = c->radius()*c->radius();
		int xmax = clat[0] + _irx; int xmin = clat[0] - _irx;
		int ymax = clat[1] + _iry; int ymin = clat[1] - _iry;
		int z_b = clat[2] - _irz;
		int z_u = clat[2] + _irz;
		int zmin = z_b >= 1 ? z_b : 1;
		int zmax = z_u < NZ ? z_u : NZ - 1;
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

			for (l = ymin; l <= ymax; l++) {
				//imy = l;
				my = l * dy;
				ipy = l;
				if (l < 0) {
					ipy += NY;
				}
				else if (l >= NY) {
					ipy -= NY;
				}

				for (m = zmin; m <= zmax; m++) {
					//ipz = imz = iz + m;
					ipz = m;
					mz = m * dz;
					diffx = mx , cv[0]();
					diffy = my - cv[1]();
					diffz = mz - cv[2]();

					if (fabs(diffx) >= 0.5*LX) {
						if (diffx > 0) {
							diffx -= LX;
						}
						else {
							diffx += LX;
						}
					}

					if (fabs(diffy) >= 0.5*LY) {
						if (diffy > 0) {
							diffy -= LY;
						}
						else {
							diffy += LY;
						}
					}

					distSq = diffx*diffx + diffy*diffy + diffz*diffz;
					if (distSq < cradSq) {
						cell_map2[ipx][ipy][ipz] = 1;
						if (distSq < normal_radSq) {
							cell_map[ipx][ipy][ipz] = c.get();
							if (c->state() == AIR) {
								air_stim_flg[ipx][ipy][ipz] = 1;
							}

							if (c->state() == ALIVE&&c->state() == DEAD) {
								age_map[ipx][ipy][ipz] = c->agek(); //‚±‚±‚Í“¯Žž‚É‘‚«‚±‚Þ‚ÆŠëŒ¯ use CAS
								cornif_map[ipx][ipy][ipz] = 1;
							}
							
						}
					}
					//field->cell_map2[ipx][ipy][ipz] = is_within_cell(mx, my, mz, cx, cy, cz, crad);
				}
			}
		}

	});
}

void Field::calc_b() {

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
		sscanf_s(line.c_str(), "%*d %d %lf %lf %lf %*f %lf %lf %lf %*f %d %lf %lf %d %lf %d %d",
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
			agek, ageb,
			ex_fat, fat,//ex_fat in_fat
			spr_len,
			0,//spring force (unused)
			state == FIX ? agki_max_fix :
			state == MUSUME ?
			agki_max : 0,//div thresh
			0,//poisson??
			div_times,//rest div times
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
	int j, jj, kk;
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