#include "field.h"
#include "cell.h"
#include <algorithm>
#include <cassert>
#include <string>
#include <iostream>
#include <sys/utime.h>



void Field::init_unused_id(Arr1D<uint32_t> &id_arr, int &next_id_var) {
	id_arr = Arr1D<uint32_t>(MAX_CELL_NUM, 0);
	next_id_var = 0;
	for (size_t i = 0; i < id_arr.size(); ++i) {
		id_arr[i] = i + 1; //1～MAX_CELL_NUMを未使用IDとする
	}
}


void Field::recycle_unused_id() {
	std::vector<bool> used_flag(MAX_CELL_NUM + 1, false); //+1は0の分
	//インデックスとidが対応
	for (auto &cellset : a_cell) {
		for (auto &c : cellset) {
			assert(c->id != 0);//0は使ってないはず
			assert(used_flag[c->id] == false);//trueならidが重複して存在しているのでおかしい
			used_flag[c->id] = true;
		}
	}
	unused_id.clear();
	for (size_t i = 1; i < used_flag.size(); ++i) {
		if (!used_flag[i]) unused_id.push_back(i);//未使用なら追加
	}
	next_id_idx = 0;
}

int Field::generate_new_id() {
	uint32_t new_id = unused_id[next_id_idx];
	if (next_id_idx == unused_id.size() - 1) {
		recycle_unused_id();
	}
	else {
		next_id_idx++;
	}
	return new_id;
}

void Field::add_future_cell(const CELL_STATE& cstate, const CellData &cdata) {
	uint32_t new_id = generate_new_id();
	auto cptr = std::shared_ptr<Cell>(new Cell(new_id, cdata));
	pending_born_cell[cstate].push_back(cptr);
	cell_by_id[new_id] = cptr;
	//alloc_id_ptr(a_cell[cstate]); //push_backでアドレスが変わるかもしれないので
}

void Field::alloc_id_ptr(const Arr1D<std::shared_ptr<Cell>>& carr) {
	for (const auto &c : carr) {
		assert(c != nullptr && c->id != 0);
		cell_by_id[c->id] = c;
	}
}

void Field::realloc_id_ptr_all() {
	std::fill(cell_by_id.begin(), cell_by_id.end(), nullptr);
	for (auto &cellset : a_cell) {
		alloc_id_ptr(cellset);
	}
}

void Field::clean_up() {
	for (auto &cellset : a_cell) {
		auto it = cellset.begin();
		while (it != cellset.end()) {
			if ((*it)->pending_kill) {
				assert(cell_by_id[(*it)->id] != nullptr);
				cell_by_id[(*it)->id] = nullptr;
				it = cellset.erase(it);
			}
			else {
				++it;
			}
		}
	}
}

void Field::born_update() {
	for (size_t i = 0; i < pending_born_cell.size(); ++i) {
		for (auto& new_cell : pending_born_cell[i]) {
			a_cell[i].push_back(new_cell);
			if (new_cell->data.div_pair_cell != nullptr) {
				new_cell->data.div_pair_cell->data.div_pair_cell = new_cell;
				paired_cell_lower_id.push_back(new_cell);
			}
			new_cell->pending_born = false;
		}
		pending_born_cell[i].clear();
	}
}

void Field::state_update() {
	for (auto &cellset : a_cell) {
		auto it = cellset.begin();
		while (it != cellset.end()) {
			if ((*it)->state_changed) {

				assert(&(a_cell[(*it)->state_dest]) != &cellset);
				a_cell[(*it)->state_dest].push_back(*it);
				(*it)->state_changed = false;
				it = cellset.erase(it);

			}
			else {
				++it;
			}
		}
	}
}

void Field::update_all_cells() {
	clean_up();
	state_update();
	born_update();
	
	update_data();
	
}

void Field::copy_cell_data_to_old() {
	for (auto &cellset : a_cell) {
		for (auto &c : cellset) {
			c->copy_data_to_old();
		}
	}
}

void Field::state_renew() {
	//FIX cell
	for (auto &cfix : a_cell[FIX]) {

		assert((cfix->old_data.dermis = find_dermis(cfix)) != nullptr);

		if (cfix->check_divide<FIX>(cfix)) {
			add_future_cell(MUSUME, cfix->div_cell_data);
		}
		else {
			cfix->data.ageb = cfix->old_data.ageb + DT_Cell*ageb_const<FIX>(cfix);
		}

	}

	for (auto &cmusume : a_cell[MUSUME]) {

		cmusume->old_data.dermis = find_dermis(cmusume);
		if (cmusume->old_data.dermis == nullptr && cmusume->old_data.div_pair_cell == nullptr) {
			if (SYSTEM == WHOLE) {
				cmusume->set_change_state(ALIVE);
			}
			else if(SYSTEM==BASAL) {
				cmusume->mark_as_delete();
			}
		}
		else {

			if (cmusume->check_divide<MUSUME>(cmusume)) {
				add_future_cell(MUSUME, cmusume->div_cell_data);
			}
			else {
				cmusume->data.ageb = cmusume->old_data.ageb + DT_Cell*ageb_const<MUSUME>(cmusume);
			}

		}
	}

	for (auto &cdead : a_cell[DEAD]) {
		if (cdead->old_data.agek >= ADHE_CONST && cdead->old_data.connected_count <= DISA_conn_num_thresh) {
			cdead->mark_as_delete();
		}
		else {
			cdead->data.agek = cdead->old_data.agek+ DT_Cell*agek_const<DEAD>(cdead);
		}
	}

	//same as DEAD
	for (auto &cair : a_cell[AIR]) {
		if (cair->old_data.agek >= ADHE_CONST && cair->old_data.connected_count <= DISA_conn_num_thresh) {
			cair->mark_as_delete();
		}
		else {
			cair->data.agek = cair->old_data.agek + DT_Cell*agek_const<AIR>(cair);
		}
	}

	for (auto &calive : a_cell[ALIVE]) {
		if (calive->old_data.agek >= THRESH_DEAD) {
			calive->set_change_state(DEAD);
			dead_count++;
		}
		else {
			calive->data.agek = calive->old_data.agek + DT_Cell*agek_const<ALIVE>(calive);
			double tmp = k_lipid_release(calive->old_data.ca2p_avg, calive->old_data.agek)
				*calive->old_data.in_fat;
			calive->data.in_fat = calive->old_data.in_fat + DT_Cell*(
				k_lipid(calive->old_data.ca2p_avg, calive->old_data.agek)*(1 - calive->old_data.in_fat)
				- tmp);
			calive->data.re_ex_fat = calive->old_data.re_ex_fat + DT_Cell*tmp;

		}
	}
}

void Field::cell_dyn() {
	//add bend
	memb_bending();
	
	
	//do not double interact


	cell_interact<MEMB,
		true, memb_memb_coef
	>();
	
	cell_interact<DER,
		true, common_common_coef,	//der_memb
		true, der_der_coef,			//der_der
		true, common_common_coef,	//der_alive
		true, common_common_coef,	//der_fix
		true, common_common_coef,	//der_musume
		true, common_common_coef,	//der_air
		true, common_common_coef	//der_dead
	>();

	//DISA do nothing
	//cell_interact<DISA>();

	cell_interact<ALIVE,
		true, supra_others_coef,	//alive-memb
		false, supra_others_coef,//alive-der NOTE:this can be ignored due to PROCESS ORDER
			//If set this true,works order-independently but runtime performance will drop
		true, supra_others_coef, //alive-alive
		true, supra_others_coef, //alive-fix
		true, supra_others_coef, //alive-musume
		true, supra_others_coef, //alive-air
		true, supra_others_coef //alive-dead
	>();

	//same as ALIVE
	cell_interact<DEAD,
		true, supra_others_coef,
		false, supra_others_coef,

		false, supra_others_coef,
		true, supra_others_coef,
		true, supra_others_coef,
		true, supra_others_coef,
		true, supra_others_coef
	>();

	//same as ALIVE
	cell_interact<AIR,
		true, supra_others_coef,
		false, supra_others_coef,
		false, supra_others_coef,
		true, supra_others_coef,
		true, supra_others_coef,
		true, supra_others_coef,
		false, supra_others_coef
	>();


	cell_interact<FIX,
		true, stem_memb_coef<FIX>,
		false, null_coef,				//order-dependent
		false, null_coef,				//order-dependent
		true, stem_stem_coef,
		true, stem_stem_coef,
		false, null_coef,
		false, null_coef
	>();

	cell_interact<MUSUME,
		true, stem_memb_coef<MUSUME>,
		false, null_coef,				//order-dependent
		false, null_coef,				//order-dependent
		false, stem_stem_coef,
		true, stem_stem_coef,
		false, null_coef,
		false, null_coef
	>();

	pair_interact();


	
	if (KEY_DERMAL_CHANGE == 0) {
		for (auto& cmemb : a_cell[MEMB]) {
			cmemb->data.x = cmemb->old_data.x;
			cmemb->data.y = cmemb->old_data.y;
			cmemb->data.z = cmemb->old_data.z;
		}
	}

	/*
		
	*/
	//cell_periodic_pos_fix();
	for (auto& cset : a_cell) {
		for (auto& c : cset) {
			c->old_data.x = c->data.x;
			c->old_data.y = c->data.y;
			c->old_data.z = c->data.z;
		}
	}
	
	state_renew();

	pair_div();

	cell_periodic_pos_fix();

	connect_cell();
}

void Field::pair_interact() {
	//pair interact
	for (auto& paired : paired_cell_lower_id) {
		assert(paired->old_data.div_pair_cell != nullptr);
		CellPtr& oppo = paired->old_data.div_pair_cell;
		assert(oppo->old_data.div_pair_cell == paired);
		cell_intr_tmpl<pair_coef>(paired, oppo);
	}
}

void Field::pair_div() {
	for (auto& paired : paired_cell_lower_id) {
		assert(paired->old_data.div_pair_cell != nullptr);
		CellPtr& oppo = paired->old_data.div_pair_cell;
		assert(oppo->old_data.div_pair_cell == paired);

		if (paired->old_data.pair_spr_nat_len < 2.0*paired->old_data.radius) {

			oppo->data.pair_spr_nat_len =
				paired->data.pair_spr_nat_len = paired->old_data.pair_spr_nat_len + DT_Cell*eps_L;

		}
		else {
			double distSq = dist3Sq(paired->old_data.x, paired->old_data.y, paired->old_data.z,
				oppo->old_data.x, oppo->old_data.y, oppo->old_data.z);
			double rad_sum = paired->old_data.radius + oppo->old_data.radius;
			double unpair_distSq = unpair_dist_coef*unpair_dist_coef*rad_sum*rad_sum;
			if (distSq > unpair_distSq) {
				oppo->data.pair_spr_nat_len =
					paired->data.pair_spr_nat_len = 0;
				oppo->data.div_pair_cell =
					paired->data.div_pair_cell = nullptr;
				oppo->data.spring_force =
					paired->data.spring_force = 0;
				paired->mark_as_got_unpaired();
			}
		}
	}
}

void Field::cell_periodic_pos_fix() {
	for (auto& cset : a_cell) {
		for (auto& c : cset) {
			if (c->data.x > LX) {
				c->data.x -= LX;
			}
			else if (c->data.x < 0) {
				c->data.x += LX;
			}
			

			if (c->data.y > LY) {
				c->data.y -= LY;
			}
			else if (c->data.y < 0) {
				c->data.y += LY;
			}
		}
	}
}

void Field::update_cell_grid() {

	for (auto& cx : cell_grid) {
		for (auto& cy : cx) {
			for (auto& cz : cy) {
				for (auto& cset : cz) {
					cset.clear();
				}
			}
		}
	}

	for (auto& cx : cell_grid_count) {
		for (auto& cy : cx) {
			for (auto& cz : cy) {
				cz = 0;
			}
		}
	}

	int aix, aiy, aiz;
	for (size_t i = 0; i < a_cell.size(); ++i) {
		if (i == DISA)continue;
		for (auto& cell : a_cell[i]) {
			aix = (int)(ret_x(LX / 2., cell->old_data.x) / AREA_GRID);
			aiy = (int)(ret_y(LY / 2., cell->old_data.y) / AREA_GRID);
			aiz= (int)((cell->old_data.z < 0. ? 0 : cell->old_data.z) / AREA_GRID);

			assert(!(aix >= ANX || aiy >= ANY || aiz >= ANZ || aix < 0 || aiy < 0 || aiz < 0));

			cell_grid[aix][aiy][aiz][i].push_back(cell);
			cell_grid_count[aix][aiy][aiz]++;

			assert(cell_grid_count[aix][aiy][aiz] < N3);
		}
	}
}

void Field::connect_cell_grid_ranged(int th_num,std::atomic<bool>* need_calc,
	CELL_STATE* cst,double* cx,double*cy,double*cz,double*crad,const CellPtr* read_only_cell,
	Arr2D<CellPtr>* out_conn_cell,int* out_conn_count,bool* out_touch) {
	int p, jj, kk, jl, jr, jb, ju;
	int ii=2;
	int aix, aiy, aiz;
	int anx, any, anz;
	double rad_sum;
	double diffx, diffy, diffz;
	bool dead_connect;
	while (1) {
		while (!th_need_calc[th_num]) {
			
			//std::this_thread::sleep_for(std::chrono::milliseconds(1));
		}
		anx = (int)(*cx / AREA_GRID);
		any = (int)(*cy / AREA_GRID);
		anz = (int)(*cz / AREA_GRID);
		dead_connect = false;
		int j = anx - ii + th_num;
		aix = (j + ANX) % ANX;
		for (int k = any - ii; k <= any + ii; k++) {
			aiy = (k + ANY) % ANY;
			for (int l = anz - ii; l <= anz + ii; l++) {
				aiz = (l + ANZ) % ANZ;
				//int csize = cell_grid_count[aix][aiy][aiz]; //avoid slow debug
				for (int state = 0; state < STATE_NUM; ++state) {
					if (state == DISA || (*cst == MEMB && state == MEMB))continue;
					for (auto& cg : cell_grid[aix][aiy][aiz][state]) {
						if (*read_only_cell == cg) continue;
						diffx = perio_diff_x(*cx, cg->old_data.x);
						diffy = perio_diff_y(*cy, cg->old_data.y);
						diffz = *cz - cg->old_data.z;
						rad_sum = *crad + cg->old_data.radius;
						if (diffx*diffx + diffy*diffy + diffz*diffz <= SQ(LJ_THRESH*rad_sum)) {
							(*out_conn_cell)[state].push_back(cg);
							(*out_conn_count)++;
							assert(*out_conn_count < N2);
							dead_connect = dead_connect || state == DEAD;
						}
					}

				}
			}
		}

		*out_touch = (*cst == ALIVE && dead_connect);
		th_need_calc[th_num] = false;
		//printf("%d\n", my_th_num);
	}
		
	
}

void Field::connect_cell() {
	assert(memb_indices_are_set);
	update_cell_grid();
	//static bool memb_connected = false;
//	int p, ii, jj, kk, jl, jr, jb, ju;
	int ii;
	int aix, aiy, aiz;
	int anx, any, anz;
	for (auto& cset : a_cell) {
		for (auto& cell : cset) {
			cell->data.connected_count = 0;
			for (auto& conn_cell_set : cell->data.connected_cell) {
				conn_cell_set.clear();
			}
		}
	}
	//need memb indices setup

	for (size_t j = 0; j < a_cell[MEMB].size(); j++) {

		auto& conn_cell = a_cell[MEMB][j]->data.connected_cell[MEMB];
		conn_cell.push_back(a_cell[MEMB][memb_l[j]]);
		conn_cell.push_back(a_cell[MEMB][memb_r[j]]);
		conn_cell.push_back(a_cell[MEMB][memb_b[j]]);
		conn_cell.push_back(a_cell[MEMB][memb_u[j]]);
		a_cell[MEMB][j]->data.connected_count += 4;
	}
	double rad_sum;
	double cx, cy, cz, crad;
	double diffx, diffy, diffz;
	for (int i = 0; i < a_cell.size(); ++i) {
		if (i == DISA)continue;
		th_cst = (CELL_STATE)i;
		for (auto& cell : a_cell[i]) {
			th_cx = cell->old_data.x;
			th_cy = cell->old_data.y;
			th_cz = cell->old_data.z;
			th_crad = cell->old_data.radius;
			th_RO_cell = cell;
			for (auto& cc : th_out_conn_cell) {
				for (auto& cset : cc) {
			
					cset.clear();
				}
			}
				for (auto& cnt : th_out_conn_count) {
					cnt = 0;
				}
			

			for (auto& t : th_out_touch) {
				t = false;
			}
			th_need_calc[0] = true;
			th_need_calc[1] = true;
			th_need_calc[2] = true;
			th_need_calc[3] = true;
			th_need_calc[4] = true;
			
			while (th_need_calc[0] || th_need_calc[1] || th_need_calc[2] || th_need_calc[3] || th_need_calc[4]) {
				//printf("aa");
			}
			
			for (int k = 0; k <= 4; k++) {
				for (int stt = 0; stt < STATE_NUM; stt++) {
					auto& v1 = cell->data.connected_cell[stt];
					auto& v2 = th_out_conn_cell[k][stt];
					
					//assert(v2.size() == 0);
					v1.insert(v1.end(), v2.begin(), v2.end());
				}
				cell->data.connected_count += th_out_conn_count[k];
				cell->data.touch = cell->data.touch || th_out_touch[k];
			}
			/*
			if (i==9 && cell->data.connected_cell[9].size() != 0) {
				printf("ueyso %d\n", cell->data.connected_cell[9].size());
			}
			*/
		}
	}
	/*
	anx = (int)(cell->old_data.x / AREA_GRID);
	any = (int)(cell->old_data.y / AREA_GRID);
	anz = (int)(cell->old_data.z / AREA_GRID);

	assert(!(anx >= ANX || any >= ANY || anz >= ANZ || anx < 0 || any < 0 || anz < 0));

	ii = 2;
	bool dead_connect = false;
	for (int j = anx - ii; j <= anx + ii; j++) {
		aix = (j + ANX) % ANX;
		for (int k = any - ii; k <= any + ii; k++) {
			aiy = (k + ANY) % ANY;
			for (int l = anz - ii; l <= anz + ii; l++) {
				aiz = (l + ANZ) % ANZ;
				//int csize = cell_grid_count[aix][aiy][aiz]; //avoid slow debug
				for (int state = 0; state < STATE_NUM; ++state) {
					if (state == DISA || (i == MEMB && state == MEMB))continue;
					for (auto& cg : cell_grid[aix][aiy][aiz][state]) {
						if (cell == cg) continue;
						diffx = perio_diff_x(cx , cg->old_data.x);
						diffy = perio_diff_y(cy , cg->old_data.y);
						diffz = cz - cg->old_data.z;
						rad_sum =crad  + cg->old_data.radius;
						if (diffx*diffx+diffy*diffy+diffz*diffz <= SQ(LJ_THRESH*rad_sum)) {
						cell->data.connected_cell[state].push_back(cg);
							cell->data.connected_count++;
							assert(cell->data.connected_count < N2);
							dead_connect = dead_connect || state == DEAD;
						}
					}

				}
			}
		}
	}
	cell->data.touch =( i == ALIVE && dead_connect);
	*/


}

void Field::update_data() {
	copy_cell_data_to_old();
}

void Field::set_memb_indices() {
	int jj, kk;
	int memb_num = (int)a_cell[MEMB].size();
	memb_u.clear(); memb_u.resize(memb_num);
	memb_b.clear(); memb_b.resize(memb_num);
	memb_l.clear(); memb_l.resize(memb_num);
	memb_r.clear(); memb_r.resize(memb_num);
	memb_ll.clear(); memb_ll.resize(memb_num);
	memb_bb.clear(); memb_bb.resize(memb_num);

	for (int j = 0; j < memb_num; j++) {
		jj = j % NMX;
		kk = j / NMX;
		if (jj == 0) memb_l[j] = j + NMX - 1;
		else memb_l[j] = j - 1;

		if (jj <= 1) memb_ll[j] = j + NMX - 2;
		else memb_ll[j] = j - 2;

		if (jj == NMX - 1) memb_r[j] = j - (NMX - 1);
		else memb_r[j] = j + 1;

		if (kk == 0) memb_b[j] = j + NMX*NMY - NMX;
		else memb_b[j] = j - NMX;

		if (kk <= 1) memb_bb[j] = j + NMX*NMY - 2 * NMX;
		else memb_bb[j] = j - 2 * NMX;

		if (kk == NMY - 1) memb_u[j] = j - (NMX*NMY - NMX);
		else memb_u[j] = j + NMX;

		if (memb_l[j] < 0 || memb_l[j] >= memb_num) { printf("error: illegal index\n"); exit(1); }
		if (memb_ll[j] < 0 || memb_ll[j] >= memb_num) { printf("error: illegal index\n"); exit(1); }
		if (memb_r[j] < 0 || memb_r[j] >= memb_num) { printf("error: illegal index\n"); exit(1); }
		if (memb_b[j] < 0 || memb_b[j] >= memb_num) { printf("error: illegal index\n"); exit(1); }
		if (memb_bb[j] < 0 || memb_bb[j] >= memb_num) { printf("error: illegal index\n"); exit(1); }
		if (memb_u[j] < 0 || memb_u[j] >= memb_num) { printf("error: illegal index\n"); exit(1); }
	}
	memb_indices_are_set = true;
}

void Field::memb_bending() {
	assert(memb_indices_are_set);
	auto& mset = a_cell[MEMB];
	int memb_num = (int)mset.size();
	static Arr1D<double> distr(memb_num),rx(memb_num), ry(memb_num), rz(memb_num);
	static Arr1D<double> distu(memb_num),ux(memb_num), uy(memb_num), uz(memb_num);
	static Arr1D<double> ipu(memb_num), ipr(memb_num);
	for (int j = 0; j < memb_num; j++) {
		auto& old = mset[j]->old_data;
		auto& oldr = mset[memb_r[j]]->old_data;
		distr[j] =
			sqrt(dist3Sq(
				oldr.x, oldr.y, oldr.z, old.x, old.y, old.z
				));
		assert(distr[j] >= 1e-5);
		rx[j] = perio_diff_x(oldr.x, old.x) / distr[j];
		ry[j] = perio_diff_y(oldr.y, old.x) / distr[j];
		rz[j] = (oldr.z - old.z) / distr[j];
		auto& oldu = mset[memb_u[j]]->old_data;
		distu[j] =
			sqrt(dist3Sq(
				oldu.x, oldu.y, oldu.z, old.x, old.y, old.z
				));
		assert(distu[j] >= 1e-5);
		ux[j] = perio_diff_x(oldu.x, old.x) / distu[j];
		uy[j] = perio_diff_y(oldu.y, old.x) / distu[j];
		uz[j] = (oldu.z - old.z) / distu[j];
	}

	for (int j = 0; j < memb_num; j++) {
		int rj = memb_r[j], uj = memb_u[j];
		ipr[j] = rx[rj] * rx[j] + ry[rj] * ry[j] + rz[rj] * rz[j];
		ipu[j] = ux[uj] * ux[j] + uy[uj] * uy[j] + uz[uj] * uz[j];
	}

	for (int j = 0; j < memb_num; j++) {
		int jr = memb_r[j], ju = memb_u[j],jl=memb_l[j],jb=memb_b[j],
			jll=memb_ll[j],jbb=memb_bb[j];

		mset[j]->data.x = mset[j]->old_data.x
			+ DT_Cell*bend_force_sqr(rx, ux, ipr, ipu, distr, distu, j, jr, jl, jll, ju, jb, jbb);
		mset[j]->data.y = mset[j]->old_data.y
			+ DT_Cell*bend_force_sqr(ry, uy, ipr, ipu, distr, distu, j, jr, jl, jll, ju, jb, jbb);
		mset[j]->data.z = mset[j]->old_data.z
			+ DT_Cell*bend_force_sqr(rz, uz, ipr, ipu, distr, distu, j, jr, jl, jll, ju, jb, jbb);
	}
}

//for disp
void Field::calc_energy() {

}

double Field::get_z_max() {
	double z = 0;
	for (int state = 0; state < a_cell.size(); state++) {
		if (state == DEAD || state == ALIVE || state == MUSUME || state == FIX) {
			for (auto& cell : a_cell[state]) {
				if (z < cell->old_data.z) {
					z = cell->old_data.z;
				}
			}
		}
	}
	return z;
}

CellPtr Field::add_new_cell_direct(const CELL_STATE& cstate, const CellData &cdata) {

	uint32_t new_id = generate_new_id();
	auto cptr = std::shared_ptr<Cell>(new Cell(new_id, cdata));
	a_cell[cstate].push_back(cptr);
	cell_by_id[new_id] = cptr;
	cptr->pending_born = false;
	return cptr;
}

void Field::init_with_file(std::ifstream& dstrm) {
	std::string line;
	CELL_STATE state;
	int div_times,touch,pair_cell_id,stem_orig_id;
	double rad, ageb, agek, x, y, z, fat, spr_len,ex_fat;
	int id_count = 0;
	std::vector<CELL_STATE> tmp_pair_state(MAX_CELL_NUM, UNUSED);
	std::vector<CellPtr> tmp_pair(MAX_CELL_NUM,nullptr);
	while (!dstrm.eof()) {
		std::getline(dstrm, line);
		sscanf_s(line.c_str(), "%*d %d %lf %lf %lf %*f %lf %lf %lf %*f %d %lf %lf %d %lf %d %d",
			&state,&rad,&ageb,&agek,&x,&y,&z,&div_times,&fat,&ex_fat,&touch,&spr_len,&pair_cell_id,&stem_orig_id);
		if (state == BLANK)return; //owari
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
		std::cout << "id:" << id_count << std::endl;
		CellData cdata;
		cdata.ageb = ageb;
		cdata.agek = agek;
		cdata.radius = rad;
		if (state == MUSUME) {
			cdata.div_age_thresh = agki_max;
		}
		else if (state == FIX) {
			cdata.div_age_thresh = agki_max_fix;
			
		}
		
		cdata.x = x;
		cdata.y = y;
		cdata.z = z;
		
		cdata.rest_div_times = div_times;
		cdata.re_ex_fat = fat;
		cdata.in_fat = ex_fat;
		cdata.touch = touch == 1;
		cdata.pair_spr_nat_len = spr_len;
		auto cptr = add_new_cell_direct(state, cdata);
		if (pair_cell_id > -1) {
			if (tmp_pair[id_count] != nullptr) {
				assert(tmp_pair[id_count]->old_data.div_pair_cell == nullptr);
				cptr->old_data.div_pair_cell = tmp_pair[id_count];
				tmp_pair[id_count]->old_data.div_pair_cell = cptr;
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
			}
			else {
				tmp_pair[pair_cell_id] = cptr;
				tmp_pair_state[pair_cell_id] = state;
			}
		}
		
		id_count++;
	}
}
class ca_dyn_init_current {
public:
	void operator()(CellPtr& c) {
		for (int i = 0; i < STATE_NUM; i++) {
			int pre_s = c->old_data.gj[i].size();
			int cell_s = c->old_data.connected_cell[i].size();
			c->old_data.gj[i].resize(cell_s);
			for (int j = pre_s; j < cell_s; j++) {
				c->old_data.gj[i][j] = w0;
			}
		}
		c->old_data.ca2p_avg = 0;
	}
};


void Field::ca_dyn() {
	proc_all_state_cell<ca_dyn_init_current,ALIVE,DEAD,FIX,MUSUME,AIR>();


	for (auto& cs : a_cell) {
		
	}
	
	for (int itr = 0; itr < Ca_ITR; itr++) {

	}
}