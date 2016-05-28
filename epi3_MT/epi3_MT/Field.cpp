#define _CRT_SECURE_NO_WARNINGS

#include "Field.h"
#include <string>
#include <cassert>
#include <iostream>
#include <istream>
#include <fstream>
#include <tbb/task.h>
#include <tbb/task_group.h>
#include <chrono>
#include <sstream>
#include <iomanip>
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
//maybe can be multithreaded
void Field::cell_state_renew() {
    cells.other_foreach_parallel([&](CellPtr& c, int i) {
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
            if(c->agek()>=cont::THRESH_DEAD){

                //sw++;
                printf("sw updated:%d\n",++sw);
            }
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
    static std::atomic<int> aindx[ANX][ANY][ANZ]={};//RawArr3D<std::atomic<int>,ANX,ANY,ANZ> aindx = { };
    static Cell* area[ANX][ANY][ANZ][N3]={nullptr};//RawArr3D<RawArr1D<Cell*, N3>, ANX, ANY, ANZ> area = { nullptr };
    //static bool mflg = false;
	//warn
	//impl-dependent-memset
	std::memset(aindx, 0, sizeof(std::atomic<int>)*ANX*ANY*ANZ);
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
		c->connected_cell.set_count(c->state() == MEMB ? 4 : 0);
	});
	cells.non_memb_foreach_parallel_native([&](CellPtr& c) {
			
			int anx = (int)(c->pos[0]() / AREA_GRID);
			int any = (int)(c->pos[1]() / AREA_GRID);
			int anz = (int)(c->pos[2]() / AREA_GRID);

			assert(!(anx >= ANX || any >= ANY || anz >= ANZ || anx < 0 || any < 0 || anz < 0));
            //Vec3<double> diffv;
			double rad_sum; double diffx, diffy, diffz;
            const int ii =2;int aix, aiy, aiz;
            int yidx[2*ii+1],zidx[2*ii+1];int yc=0,zc=0;
            for (int k = any - ii; k <= any + ii; k++,yc++){
                aiy = k;
                if (k < 0) {
                    aiy += ANY;
                }
                else if (k >= ANY) {
                    aiy -= ANY;
                }
                yidx[yc]=aiy;
            }
            for (int l = anz - ii; l <= anz + ii; l++,zc++) {
                aiz= l;
                if (l < 0) {
                    aiz += ANZ;
                }
                else if (l >= ANZ) {
                    aiz -= ANZ;
                }
                zidx[zc]=aiz;
            }
			for (int j = anx - ii; j <= anx + ii; j++) {
				aix = j;
				if (j < 0) {
					aix += ANX;
				}
				else if (aix >= ANX) {
					aix -= ANX;
				}

                for (int k = 0; k <2*ii+1; k++) {
                    aiy=yidx[k];
                    for (int l = 0;l<2*ii+1; l++) {
                        aiz=zidx[l];
                        int sz=aindx[aix][aiy][aiz];
                        for (int m = 0; m < sz; ++m) {
							Cell* o = area[aix][aiy][aiz][m];
                            if (c->my_construction_count <= o->my_construction_count)continue;

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
								o->connected_cell.add(c);
								/*
								if (o->state() == MEMB) {
									printf("eusyo");
								}
								*/
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
        if (get_state_mask(st)&( DEAD_M | ALIVE_M | MUSUME_M |FIX_M)) {
			if (zmax < c->pos[2]())zmax = c->pos[2]();
		}
	});
	return zmax;
}

void Field::cell_dynamics() {
    tbb::task_group t;
    t.run([&]{


    cells.memb_foreach_parallel_native([](CellPtr& memb) {
		memb->memb_bend_calc1();
	});

    cells.memb_foreach_parallel_native([](CellPtr& memb) {
		memb->memb_bend_calc2();
	});

    cells.memb_foreach_parallel_native([](CellPtr& memb) {
		memb->memb_bend_interact();
	});
 });
    t.run([&]{
        interact_cell();
    });
t.wait();
cells.foreach_parallel([](CellPtr& c,int i) {
	c->pos[0].update();
	c->pos[1].update();
	c->pos[2].update();
});

	cell_state_renew();
	cells.update();
	cells.all_cell_update();

    cells.other_foreach_parallel_native([](CellPtr& c) {
        if (c->pair != nullptr&&c->my_construction_count>c->pair->my_construction_count) {
			c->pair_disperse();

		}
	});

	cell_pos_periodic_fix();
	cells.foreach_parallel([](CellPtr& c,int i) {
		c->pos[0].update();
		c->pos[1].update();
		c->pos[2].update();
		c->spring_nat_len.update();
	});
	connect_cells();
}

void Field::initialize_sc(){
    cells.other_foreach([&](CellPtr& c,int i){
if(c->state()==ALIVE&&zzmax-c->pos[2]()<8*cont::R_max){
c->agek.force_set_next_value(c->agek()>cont::THRESH_DEAD?c->agek():cont::THRESH_DEAD);
c->state=AIR;
c->agek.update();
c->state.update();
}
    });
}
void Field::check_localization(){
    cells.foreach([&](CellPtr& c,int i){
       bool touch=false;
       if(c->state()==ALIVE){
           for(Cell* conn:c->connected_cell._cell()){
               if(conn->state()==DEAD){
                   touch=true;
                   break;
               }
           }
       }
	   c->is_touch = touch;

    });
}

void Field::output_data_impl(const std::string& filename){

    std::ofstream wfile;
    wfile.open(filename,std::ios::out);
    cells.foreach([](CellPtr& c,int i){
    c->__my_index=i;
    });
    using namespace std;
    cells.foreach([&](CellPtr& c,int i){
        wfile<<i<<" "
            <<c->state()<<" "
           <<fixed<<setprecision(15)<<c->radius()<<" "
          <<fixed<<setprecision(15)<<c->ageb()<<" "
         <<fixed<<setprecision(15)<<c->agek()<<" "
            <<fixed<<setprecision(15)<<c->ca2p()<<" "
           <<fixed<<setprecision(15)<<c->pos[0]()<<" "
          <<fixed<<setprecision(15)<<c->pos[1]()<<" "
         <<fixed<<setprecision(15)<<c->pos[2]()<<" "
        <<fixed<<setprecision(15)<<c->ca2p_avg()<<" "
           <<c->rest_div_times()<<" "
          <<fixed<<setprecision(15)<<c->ex_fat()<<" "
         <<fixed<<setprecision(15)<<c->in_fat()<<" "
        <<(c->is_touch?1:0)<<" "
           <<fixed<<setprecision(15)<<c->spring_nat_len()<<" "
          <<(c->pair==nullptr?-1:c->pair->__my_index)<<" "
                                                  <<0<<std::endl;
    });
}

void Field::output_data(int idx){
std::ostringstream _fname;
_fname<< output_dir<<"/"<<idx;
output_data_impl(_fname.str());
//std::string filename=_fname.str();
/*
std::ofstream wfile;
wfile.open(filename,std::ios::out);
cells.foreach([](CellPtr& c,int i){
c->__my_index=i;
});
cells.foreach([&](CellPtr& c,int i){
    wfile<<i<<" "
        <<c->state()<<" "
       <<c->radius()<<" "
      <<c->ageb()<<" "
     <<c->agek()<<" "
        <<c->ca2p()<<" "
       <<c->pos[0]()<<" "
      <<c->pos[1]()<<" "
     <<c->pos[2]()<<" "
    <<c->ca2p_avg()<<" "
       <<c->rest_div_times()<<" "
      <<c->ex_fat()<<" "
     <<c->in_fat()<<" "
    <<(c->is_touch?1:0)<<" "
       <<c->spring_nat_len()<<" "
      <<(c->pair==nullptr?-1:c->pair->__my_index)<<" "
                                              <<0<<std::endl;
});
*/
}

void Field::main_loop()
{
    if(flg_forced_sc)printf("sc forced.\n");
    auto start=std::chrono::system_clock::now();
	for (int i = 0; i < cont::NUM_ITR; i++) {
        if(i%100==0){
            auto dur=std::chrono::system_clock::now()-start;
            printf("loop:%d elapsed[sec]:%lf\n", i,0.001*std::chrono::duration_cast<std::chrono::milliseconds>(dur).count());
        }
        if(i%CUT==0){

            printf("logging...\n");
            output_data(i/CUT);
            printf("done.\n");
        if(i%(CUT*10)==0){
            printf("saving field data...\n");
            output_field_data(i);
            printf("done.\n");
        }
        }

		cell_dynamics();
        tbb::task_group t;
        t.run([&]{
            zzmax = calc_zzmax();
        });
        t.run([&]{
            setup_map();
        });
        t.wait();
        //set_cell_lattice();
        //setup_map(); //lattice set
		calc_b();
		//b_update();

        if(i*cont::DT_Cell>cont::T_TURNOVER&&flg_forced_sc){
flg_forced_sc=false;
printf("forced cornif\n");
initialize_sc();
num_sc=cont::NUM_SC_INIT;
        }
        if (sw>=cont::SW_THRESH || num_sc>0) {
			printf("calc ca...\n");
            calc_ca();
            //cells.all_cell_update();
			printf("end.\n");
            if(num_sc>0)num_sc--;
            sw=0;
		}
	}
}

void Field::setup_map()
{
	using namespace cont;
    tbb::parallel_for(tbb::blocked_range<int>(0, NX+1), [&](const tbb::blocked_range< int >& range) {
		std::memset(&(cell_map[range.begin()]), 0, sizeof(Cell*)*range.size()*(NY + 1)*(NZ + 1));
		std::memset(&(cell_map2[range.begin()]), 0, sizeof(uint_fast8_t)*range.size()*(NY + 1)*(NZ + 1));
    });
	cells.foreach_parallel_native([&](CellPtr& c) {
		
		auto& cv = c->pos;
        c->set_lattice();
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

				
				for (zc=0,ipz = zmin; ipz <= zmax; ipz++,zc++) {
					//ipz = imz = iz + m;
					//ipz = m;
					if ((distSq = diffxSq + a_diffySq[yc] + a_diffzSq[zc]) < cradSq) {
                        if(get_state_mask(c->state())&(ALIVE_M|FIX_M|MUSUME_M)){
						cell_map2[ipx][ipy][ipz] = 1;
                }
						
						if (distSq < normal_radSq) {
							cell_map[ipx][ipy][ipz] = c.get();
						}
						
					}
				}
			}
		}

	});

    for(int l=0;l<NZ;l++){
        for(int j=0;j<NX;j++){
            cell_map2[j][NY][l]=cell_map2[j][0][l];
            cell_map[j][NY][l]=cell_map[j][0][l];
        }

        for(int k=0;k<=NY;k++){
            cell_map2[NX][k][l]=cell_map2[0][k][l];
            cell_map[NX][k][l]=cell_map[0][k][l];
        }
    }
}

void Field::calc_b() {
	using namespace cont;
	int iz_bound = (int)((zzmax + FAC_MAP*R_max) / dz);
	int* a_prev_z = new int[iz_bound];
	int* a_next_z = new int[iz_bound];
	std::swap(old_ext_stim, _ext_stim);
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
        a_prev_z[l]=prev_z;
        a_next_z[l]=next_z;
    }
	tbb::parallel_for(tbb::blocked_range<int>(0, NX ), [&](const tbb::blocked_range< int >& range) {

		for (int j = range.begin(); j != range.end(); ++j) {
            int prev_x = per_x_prev_idx[j];
            int next_x = per_x_next_idx[j];
			for (int k = 0; k < NY; k++) {
                int prev_y = per_y_prev_idx[k];
                int next_y = per_y_next_idx[k];
				for (int l = 0; l < iz_bound; l++) {
                    int prev_z = a_prev_z[l], next_z = a_next_z[l];
					double dum_age = 0;
					bool flg_cornified = false;
					
					if (cell_map[j][k][l]!=nullptr) {
						if (get_state_mask(cell_map[j][k][l] ->state())&(ALIVE_M|DEAD_M)) {
							dum_age = cell_map[j][k][l]->agek();
							flg_cornified = true;
						}
					}
					auto& cext = (*old_ext_stim)[j][k][l];
                    (*_ext_stim)[j][k][l]= (*old_ext_stim)[j][k][l]
                        +DT_Ca*(DB *
                         (cell_map2[prev_x][k][l] * ((*old_ext_stim)[prev_x][k][l] - cext)
						+ cell_map2[j][prev_y][l] * ((*old_ext_stim)[j][prev_y][l] - cext)
						+ cell_map2[j][k][prev_z] * ((*old_ext_stim)[j][k][prev_z] - cext)
						+ cell_map2[next_x][k][l] * ((*old_ext_stim)[next_x][k][l] - cext)
						+ cell_map2[j][next_y][l] * ((*old_ext_stim)[j][next_y][l] - cext)
						+ cell_map2[j][k][next_z] * ((*old_ext_stim)[j][k][next_z] - cext)) *inv_dz*inv_dz
						+ fB(dum_age, (*old_ext_stim)[j][k][l] ,flg_cornified));


				}
			}
		}
	});
	for (int l = 0; l <= iz_bound; l++) {
		for (int j = 0; j < NX; j++) (*_ext_stim)[j][NY][l]= (*_ext_stim)[j][0][l];
		for (int k = 0; k <= NY; k++)(*_ext_stim)[NX][k][l]= (*_ext_stim)[0][k][l];
	}
	delete a_prev_z;
	delete a_next_z;
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

    int* a_prev_z = new int[iz_bound];
    int* a_next_z = new int[iz_bound];
    double dummy_diffu=0;
    tbb::parallel_for(tbb::blocked_range<int>(0, NX), [&](const tbb::blocked_range< int >& range) {
        for (int j = range.begin(); j!= range.end(); ++j) {
            for (int k = 0; k < NY; k++) {
                for (int l = 0; l <iz_bound; l++) {
                    double*	tmp = &dummy_diffu;
                    uint_fast8_t asf=0;
                    if (cell_map[j][k][l] != nullptr) {
                        if (get_state_mask(cell_map[j][k][l]->state())&(ALIVE_M | FIX_M | MUSUME_M)) {
                            tmp = &(cell_map[j][k][l]->diffu);
                        }
                        int c_count= cell_map[j][k][l]->connected_cell.count();
                        auto& ccell=cell_map[j][k][l]->connected_cell._cell();
                        for(int cc=0;cc<c_count;++cc){
                            if(ccell[cc]->state()==AIR){
                                asf=1;
                                break;
                            }
                        }

                    }
                    cell_diffu_map[j][k][l]=tmp;
                    air_stim_flg[j][k][l]=asf;
                    //air_stim_flg‚Í[0,NX)*[0,NY)*[0,iz_bound)‚Å‚µ‚©Žg‚í‚ê‚È‚¢
                }
            }
        }
    });

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
        a_prev_z[l]=prev_z;
        a_next_z[l]=next_z;
    }
	for (int cstp = 0; cstp <Ca_ITR; cstp++) {
   // std::memcpy(old_ATP, _ATP, sizeof(double)*(NX + 1)*(NY + 1)*(NZ + 1));
		std::swap(old_ATP, _ATP);
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
        });
        cells.other_foreach_parallel_native([&iz_bound,this](CellPtr& c) {
            if (get_state_mask(c->state())&(ALIVE_M | FIX_M | MUSUME_M)) {
				int ix = c->lat[0]; int iy = c->lat[1]; int iz = c->lat[2];

				assert(iz < iz_bound);

				double tmp_a = grid_avg8_n(*old_ATP, ix, iy, iz);
				double tmp_B = grid_avg8_n(*_ext_stim, ix, iy, iz);

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
                c->ca2p_avg += c->ca2p.get_next_value();
			}
		});

		tbb::parallel_for(tbb::blocked_range<int>(0, NX), [&](const tbb::blocked_range< int >& range) {
			for (int j = range.begin(); j!= range.end(); ++j) {
                int prev_x = per_x_prev_idx[j];
                int next_x = per_x_next_idx[j];
				for (int k = 0; k < NY; k++) {
                    int prev_y = per_y_prev_idx[k];
                    int next_y = per_y_next_idx[k];
					for (int l = 0; l < iz_bound; l++) {
                        int prev_z = a_prev_z[l], next_z = a_next_z[l];
                        double 	tmp = *cell_diffu_map[j][k][l];
                        double catp=(*old_ATP)[j][k][l];
						(*_ATP)[j][k][l]
                            =catp+ DT_Ca*(Da * (cell_map2[prev_x][k][l] * ((*old_ATP)[prev_x][k][l] - catp)
                                + cell_map2[next_x][k][l] * ((*old_ATP)[next_x][k][l] - catp)
                                + cell_map2[j][prev_y][l] * ((*old_ATP)[j][prev_y][l] - catp)
                                + cell_map2[j][next_y][l] * ((*old_ATP)[j][next_y][l] - catp)
                                + cell_map2[j][k][prev_z] * ((*old_ATP)[j][k][prev_z] - catp)
                                + cell_map2[j][k][next_z] * ((*old_ATP)[j][k][next_z] - catp)) *inv_dx*inv_dx
                                + fa(tmp, catp)+air_stim_flg[j][k][l] * AIR_STIM);
					}
				}
			}
		});
		for (int l = 0; l <= iz_bound; l++) {
			for (int j = 0; j < NX; j++) (*_ATP)[j][NY][l]= (*_ATP)[j][0][l];
			for (int k = 0; k <= NY; k++) (*_ATP)[NX][k][l] = (*_ATP)[0][k][l];
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
	}
	cells.other_foreach_parallel_native([](CellPtr& c) {
		if (get_state_mask(c->state())&(ALIVE_M | FIX_M | MUSUME_M)) {
			c->ca2p_avg /= Ca_ITR;
		}
		else {
			c->ca2p_avg.force_set_next_value(0);
			c->ca2p.force_set_next_value(0);
		}
        c->ca2p_avg.update();
        c->ca2p.update();
	});
    delete a_prev_z;delete a_next_z;
}

void Field::init_with_file(std::ifstream& dstrm,bool use_last) {
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
    while (std::getline(dstrm, line)) {

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
            state == FIX ? agki_max_fix :
			state == MUSUME ?
#ifdef AGE_DBG
			agki_max*0.01 :
#else
		agki_max:
#endif
			0,//div thresh
            state==FIX?div_max:div_times,//rest div times
			stem_orig_id<MALIG_NUM,//malignant
			touch == 1//touch
			);
		cells.add_direct(cptr);

		if (pair_cell_id > -1) {
			if (tmp_pair[id_count] != nullptr) {
				assert(tmp_pair[id_count]->pair == nullptr);
				cptr->pair = tmp_pair[id_count];
				tmp_pair[id_count]->pair = cptr;
			}
			else {
				tmp_pair[pair_cell_id] = cptr;
				tmp_pair_state[pair_cell_id] = state;
			}
		}
printf("Phase %d  Cell loaded:%d\n",phase, id_count++);

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

    printf("initializing values...\n"); //initial_u
    cells.foreach_parallel_native([&](CellPtr& c){
c->ca2p.force_set_next_value(u0);
c->ex_inert.force_set_next_value(v0);
c->IP3.force_set_next_value(p0);
c->connected_cell.foreach([&](Cell* conn){
    c->gj.at(conn).force_set_next_value(w0);
    c->gj.at(conn).update();
});

c->ca2p.update();
c->ex_inert.update();
c->IP3.update();
    });

    cells.foreach_parallel_native([&](CellPtr& c){
    c->ca2p_avg.force_set_next_value(c->ca2p());
    c->ca2p_avg.update();
    });

    cells.other_foreach_parallel_native([&](CellPtr& c){
        if(c->state() == DEAD){
        c->ca2p_avg.force_set_next_value(0);
        c->ca2p.force_set_next_value(0);
        c->ca2p_avg.update();
        c->ca2p.update();
        }
    });
    printf("done.\n");

if(use_last){
init_field_data_with_file("last_data_uvp","last_data_w_alt","last_data_a","last_data_B");
}

    printf("Ca2P iteration loop:%d\n",Ca_ITR);
}

void Field::output_field_data(int i){
    output_data_impl("last_data_cell");

    cells.foreach([](CellPtr& c,int i){
    c->__my_index=i;
    });
FILE *fuvp,*fw_alt,*fa,*fB,*finfo; //dont use ofstream
finfo=std::fopen("last_data_info.txt","w");
std::fprintf(finfo,"i=%d,t=%lf,cell_num=%d\n",i,cont::DT_Cell*i,cells._raw_cell_set().size());
std::fclose(finfo);
tbb::task_group t;
t.run([&]{
fuvp=std::fopen("last_data_uvp","w");
cells.foreach([&fuvp](CellPtr& c,int i){
std::fprintf(fuvp,"%d %1.14e %1.14e %1.14e\n",i,c->ca2p(),c->ex_inert(),c->IP3());
});
std::fclose(fuvp);
});
t.run([&]{
fw_alt=std::fopen("last_data_w_alt","w");
cells.foreach([&fw_alt](CellPtr& c,int i){
std::fprintf(fw_alt,"%d %d\n",i,c->gj.size());
for(auto& gjv:c->gj){
    std::fprintf(fw_alt,"%d %1.14e\n",gjv.first->__my_index,gjv.second());
}
});
std::fclose(fw_alt);
});

t.run([&]{
fa=std::fopen("last_data_a","w");
for(int i=0;i<=cont::NX;i++){
    for(int j=0;j<=cont::NY;j++){
        for(int k=0;k<=cont::NZ;k++){
            std::fprintf(fa,"%f %f %f %1.14e\n",cont::dx*i,cont::dy*j,cont::dz*k,(*_ATP)[i][j][k]);
        }
        //std::fprintf(fa,"\n");
    }
    //std::fprintf(fa,"\n");
}
std::fclose(fa);
});

t.run([&]{
fB=std::fopen("last_data_B","w");
for(int i=0;i<=cont::NX;i++){
    for(int j=0;j<=cont::NY;j++){
        for(int k=0;k<=cont::NZ;k++){
            std::fprintf(fB,"%f %f %f %1.14e\n",cont::dx*i,cont::dy*j,cont::dz*k,(*_ext_stim)[i][j][k]);
        }
        //std::fprintf(fB,"\n");
    }
    //std::fprintf(fB,"\n");
}
std::fclose(fB);
});
t.wait();
}

void Field::init_field_data_with_file(std::string nuvp, std::string nw_alt, std::string na, std::string nB){
std::ifstream fuvp(nuvp),fw_alt(nw_alt),fa(na),fB(nB);
std::string ln;
int cidx,lcount;
double ca2p,ex_inert,IP3;
lcount=0;
while(std::getline(fuvp,ln)){
    sscanf(ln.c_str(),"%d %lf %lf %lf",&cidx,&ca2p,&ex_inert,&IP3);

CellPtr& c = cells._raw_cell_set()[cidx];
c->ca2p.force_set_next_value(ca2p);c->ex_inert.force_set_next_value(ex_inert);c->IP3.force_set_next_value(IP3);
c->ca2p.update();c->ex_inert.update();c->IP3.update();
lcount++;
}
if(lcount!=cells._raw_cell_set().size()){
    printf("cell num mismatch. lcount=%d real=%d\n",lcount,cells._raw_cell_set().size());
    exit(1);
}

lcount=0;
int gj_len=0;
while(std::getline(fw_alt,ln)){
if(sscanf(ln.c_str(),"%d %d",&cidx,&gj_len) != 2){
    fprintf(stdout,"bad GJ data format (index)\n");
    fflush(stdout);
    exit(-1);
}
CellPtr& c = cells._raw_cell_set()[cidx];
int gj_idx=0;double gj_val=0;
for(int i=0;i<gj_len;i++){
if(!std::getline(fw_alt,ln)){
fprintf(stdout,"bad GJ data format\n");
fflush(stdout);
exit(-1);
}
sscanf(ln.c_str(),"%d %lf",&gj_idx,&gj_val);
c->gj.emplace(cells._raw_cell_set()[gj_idx].get(),gj_val);
}
lcount++;
}
assert(lcount==cells._raw_cell_set().size());
for(int i=0;i<=cont::NX;i++){
    for(int j=0;j<=cont::NY;j++){
        for(int k=0;k<=cont::NZ;k++){
            if(!std::getline(fa,ln)){
                fprintf(stdout,"bad ATP data format\n");
                fflush(stdout);
                exit(-1);
            }
            sscanf(ln.c_str(),"%*f %*f %*f %lf",&((*_ATP)[i][j][k]));
        }
    }
}

for(int i=0;i<=cont::NX;i++){
    for(int j=0;j<=cont::NY;j++){
        for(int k=0;k<=cont::NZ;k++){
            if(!std::getline(fB,ln)){
                fprintf(stdout,"bad ext_stim data format\n");
                fflush(stdout);
                exit(-1);
            }
            sscanf(ln.c_str(),"%*f %*f %*f %lf",&((*_ext_stim)[i][j][k]));
        }
    }
}
}
