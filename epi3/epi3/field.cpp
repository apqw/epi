#include "field.h"
#include "cell.h"
#include "field_func.h"
#include <algorithm>
#include <cassert>
#include <string>
#include <iostream>
#include <sys/utime.h>
#ifdef  _WIN32
#include <Windows.h>
#endif //  _WIN32

/*
	system
*/


void Field::init_with_file(std::ifstream& dstrm) {
	std::string line;
	CELL_STATE state;
	int div_times, touch, pair_cell_id, stem_orig_id;
	double rad, ageb, agek, x, y, z, fat, spr_len, ex_fat;
	int id_count = 0;
	std::vector<CELL_STATE> tmp_pair_state(MAX_CELL_NUM, UNUSED);
	std::vector<CellPtr> tmp_pair(MAX_CELL_NUM, nullptr);
	while (!dstrm.eof()) {
		std::getline(dstrm, line);
		sscanf_s(line.c_str(), "%*d %d %lf %lf %lf %*f %lf %lf %lf %*f %d %lf %lf %d %lf %d %d",
			&state, &rad, &ageb, &agek, &x, &y, &z, &div_times, &fat, &ex_fat, &touch, &spr_len, &pair_cell_id, &stem_orig_id);
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


/*
	Cell management
*/

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

CellPtr Field::add_new_cell_direct(const CELL_STATE& cstate, const CellData &cdata) {

	uint32_t new_id = generate_new_id();
	auto cptr = std::shared_ptr<Cell>(new Cell(new_id, cdata));
	a_cell[cstate].push_back(cptr);
	cell_by_id[new_id] = cptr;
	cptr->pending_born = false;
	return cptr;
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



/*
	Cell dynamics
*/

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
	//if(i*DT > T_TURNOVER && flg_forced_sc) { ...... in 21-der.cの部分

	if (c_step*DT_Cell > T_TURNOVER && is_sc_forced) {
		for (auto &calive : a_cell[ALIVE]) {
			if (zzmax - calive->old_data.z < 8 * R_max && !calive->state_changed) {
				calive->data.agek = calive->data.agek>THRESH_DEAD ? calive->data.agek : THRESH_DEAD;
				calive->set_change_state(AIR);
			}
		}
		is_sc_forced = false;
		num_sc = NUM_SC_INIT;
	}
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

void Field::memb_bending() {
	assert(memb_indices_are_set);
	auto& mset = a_cell[MEMB];
	int memb_num = (int)mset.size();
	static Arr1D<double> distr(memb_num), rx(memb_num), ry(memb_num), rz(memb_num);
	static Arr1D<double> distu(memb_num), ux(memb_num), uy(memb_num), uz(memb_num);
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
		int jr = memb_r[j], ju = memb_u[j], jl = memb_l[j], jb = memb_b[j],
			jll = memb_ll[j], jbb = memb_bb[j];

		mset[j]->data.x = mset[j]->old_data.x
			+ DT_Cell*bend_force_sqr(rx, ux, ipr, ipu, distr, distu, j, jr, jl, jll, ju, jb, jbb);
		mset[j]->data.y = mset[j]->old_data.y
			+ DT_Cell*bend_force_sqr(ry, uy, ipr, ipu, distr, distu, j, jr, jl, jll, ju, jb, jbb);
		mset[j]->data.z = mset[j]->old_data.z
			+ DT_Cell*bend_force_sqr(rz, uz, ipr, ipu, distr, distu, j, jr, jl, jll, ju, jb, jbb);
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


/*
	Cell geometry
*/


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
//	int p, jj, kk, jl, jr, jb, ju;
	int ii=2;
	int aix, aiy, aiz;
	int anx, any, anz;
	double rad_sum;
	double diffx, diffy, diffz;
	bool dead_connect;
	while (1) {
		while (!th_need_calc[th_num]) {
#ifdef _WIN32
			Sleep(0);
#endif // _WIN32

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
	int p, jj, kk, jl, jr, jb, ju;
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
		for (int i = 0; i < a_cell.size(); i++) {
			for (auto& cell : a_cell[i]) {
				cx = cell->old_data.x;
				cy = cell->old_data.y;
				cz = cell->old_data.z;
				crad = cell->old_data.radius;
				anx = (int)(cx / AREA_GRID);
				any = (int)(cy / AREA_GRID);
				anz = (int)(cz / AREA_GRID);

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
									diffx = perio_diff_x(cx, cg->old_data.x);
									diffy = perio_diff_y(cy, cg->old_data.y);
									diffz = cz - cg->old_data.z;
									rad_sum = crad + cg->old_data.radius;
									if (diffx*diffx + diffy*diffy + diffz*diffz <= SQ(LJ_THRESH*rad_sum)) {
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
				cell->data.touch = (i == ALIVE && dead_connect);
			}
		}

//	*out_touch = (*cst == ALIVE && dead_connect);
//	double rad_sum;
//	double cx, cy, cz, crad;
//	double diffx, diffy, diffz;
	/*
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

					v1.insert(v1.end(), v2.begin(), v2.end());
				}
				cell->data.connected_count += th_out_conn_count[k];
				cell->data.touch = cell->data.touch || th_out_touch[k];
			}
		}
	}
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

void Field::set_cell_lattice() {
	for (auto& cs : a_cell) {
		for (auto& c : cs) {
			c->data.latx = ((int)(c->data.x / dx) + NX) % NX;
			c->data.laty = ((int)(c->data.y / dy) + NY) % NY;
			c->data.latz = ((int)(c->data.z / dz) + NZ) % NZ;

			assert(c->data.latx >= 0 && c->data.laty >= 0 && c->data.latz >= 0);
		}
	}
}

/*
	Field map
*/

void Field::setup_map() {
	for (auto& mx : cell_map) {
		for (auto& my : mx) {
			std::fill(my.begin(), my.end(), nullptr);
		}
	}

	for (auto& mx : cell_map2) {
		for (auto& my : mx) {
			std::fill(my.begin(), my.end(), 0);
		}
	}

	for (auto& mx : air_stim) {
		for (auto& my : mx) {
			std::fill(my.begin(), my.end(), 0);
		}
	}

	proc_all_state_cell<map_setter, ALIVE, DEAD, FIX, DER, MUSUME, AIR, MEMB>(a_cell);


	for (int l = 0; l<NZ; l++) {
		for (int j = 0; j<NX; j++) {
			cell_map[j][NY][l] = cell_map[j][0][l];
			cell_map2[j][NY][l] = cell_map2[j][0][l];
		}
		for (int k = 0; k <= NY; k++) {
			cell_map[NX][k][l] = cell_map[0][k][l];
			cell_map2[NX][k][l] = cell_map2[0][k][l];
		}
	}
}

class Field::map_setter {
public:
	template<CELL_STATE state>
	void operator()(CellPtr& c, const Field* _field) {
		Field* field = const_cast<Field*>(_field);
		double cx = c->old_data.x;
		double cy = c->old_data.y;
		double cz = c->old_data.z;
		//double crad = c->old_data.radius;
		int ix = c->old_data.latx;
		int iy = c->old_data.laty;
		int iz = c->old_data.latz;
		/*
		int _irx = irx;
		int _iry = iry;
		int _irz = irz;
		*/
		int  k, l, m, ipx, ipy, ipz, imx, imy;
		double mx, my, mz;
		/*
		int xmax = ix + _irx; int xmin = ix - _irx;
		int ymax = iy + _iry; int ymin = iy - _iry;
		int z_b = iz - _irz;
		int z_u = iz + _irz;
		int zmin = z_b >= 1 ? z_b : 1;
		int zmax = z_u < NZ ? z_u : NZ - 1;
		*/
		/*
		for (k =xmin; k <= xmax; k++) {
		//imx = k;
		mx = k * dx;
		ipx = k + (k < 0 ? NX : k >= NX ? -NX : 0);

		for (l =ymin; l <= ymax; l++) {
		//imy = l;
		my = l * dy;
		ipy = l + (l < 0 ? NY : l >= NY ? -NY : 0);

		for (m = zmin; m <= zmax; m++) {
		//ipz = imz = iz + m;
		ipz = m;
		mz = m * dz;

		if (true) {
		field->cell_map[ipx][ipy][ipz] = &c;
		STATIC_IF(state == AIR) {
		field->air_stim[ipx][ipy][ipz] = AIR_STIM;
		}STATIC_ENDIF;

		STATIC_IF(state == ALIVE || state == DEAD) {
		field->age_map[ipx][ipy][ipz] = c->old_data.agek;
		field->cornif_map[ipx][ipy][ipz] = 1;
		}STATIC_ENDIF;
		}

		}
		}
		}
		*/
		int _irx = 2 * irx;
		int _iry = 2 * iry;
		int _irz = 2 * irz;
		double diffx, diffy, diffz, distSq;
		double crad = FAC_MAP*c->old_data.radius;
		double cradSq = crad*crad;
		double normal_radSq = c->old_data.radius*c->old_data.radius;
		int xmax = ix + _irx; int xmin = ix - _irx;
		int ymax = iy + _iry; int ymin = iy - _iry;
		int z_b = iz - _irz;
		int z_u = iz + _irz;
		int zmin = z_b >= 1 ? z_b : 1;
		int zmax = z_u < NZ ? z_u : NZ - 1;
		for (k = xmin; k <= xmax; k++) {
			//imx = k;
			mx = k * dx;

			ipx = k;

		
			if (k < 0) {
				ipx += NX;
			}else if (k >= NY) {
				ipx -= NX;
			}

			for (l = ymin; l <= ymax; l++) {
				//imy = l;
				my = l * dy;
				ipy = l;
				if (l < 0) {
					ipy += NY;
				}else if (l >= NY) {
					ipy -= NY;
				}

				for (m = zmin; m <= zmax; m++) {
					//ipz = imz = iz + m;
					ipz = m;
					mz = m * dz;
					diffx = mx - cx;
					diffy = my - cy;
					diffz = mz - cz;

					if (fabs(diffx) >= 0.5*LX) {
						if (diffx>0) {
							diffx -= LX;
						}
						else {
							diffx += LX;
						}
					}

					if (fabs(diffy) >= 0.5*LY) {
						if (diffy>0) {
							diffy -= LY;
						}
						else {
							diffy += LY;
						}
					}

					distSq = diffx*diffx + diffy*diffy + diffz*diffz;
					if (distSq < cradSq) {
						field->cell_map2[ipx][ipy][ipz] = 1;
						if (distSq < normal_radSq) {
							field->cell_map[ipx][ipy][ipz] = &c;
							STATIC_IF(state == AIR) {
								field->air_stim[ipx][ipy][ipz] = AIR_STIM;
							}STATIC_ENDIF;

							STATIC_IF(state == ALIVE || state == DEAD) {
								field->age_map[ipx][ipy][ipz] = c->old_data.agek;
								field->cornif_map[ipx][ipy][ipz] = 1;
							}STATIC_ENDIF;
						}
					}
					//field->cell_map2[ipx][ipy][ipz] = is_within_cell(mx, my, mz, cx, cy, cz, crad);
				}
			}
		}


	}
};


/*
	Ca (and ext stim) dynamics
*/

class Field::ca_dyn_init_current {
public:
	template<CELL_STATE state>
	void operator()(CellPtr& c, const Field* field) {
		for (int i = 0; i < STATE_NUM; i++) {
			int pre_s = c->old_data.gj[i].size();
			int cell_s = c->old_data.connected_cell[i].size();
			c->old_data.gj[i].resize(cell_s);
			
			for (int j = pre_s; j < cell_s; j++) {
				c->old_data.gj[i][j] = w0;
			}
			c->data.gj[i] = c->old_data.gj[i];
		
		}
		c->data.ca2p = c->old_data.ca2p;
		c->data.ex_inert = c->old_data.ex_inert;
		c->data.IP3 = c->old_data.IP3;
		c->old_data.ca2p_avg =c->data.ca2p_avg= 0;
		Field* fld = const_cast<Field*>(field);
		fld->ATP = fld->old_ATP;

	}
};

class Field::ca_dead_proc {
public:
	template<CELL_STATE state>
	void operator()(const CellPtr& c, const Field* _field) {
		//c->data.IP3 = c->old_data.IP3 -DT_Ca*Kpp*c->old_data.IP3

		c->data.IP3 = c->old_data.IP3*(1.0 - DT_Ca*Kpp);
		for (auto& cc : c->old_data.connected_cell[ALIVE]) {
			c->data.IP3 = c->old_data.IP3 + DT_Ca*dp*(cc->old_data.IP3 - c->old_data.IP3);
		}

		
	}
};

class Field::ca_common_cell_proc_conn {
private:
	CellPtr conn_cell;
public:
	ca_common_cell_proc_conn(CellPtr& cc) {
		conn_cell = cc;
	};
	template<CELL_STATE state>
	void operator()(CellPtr& c, const Field* _field) {
		conn_cell->data.ca2p_diff = conn_cell->old_data.ca2p_diff +
			ca2p_du*IAG<state>(conn_cell->old_data.agek)*conn_cell->old_data.gj[]
	}
};

class Field::ca_common_cell_proc {
public:
	template<CELL_STATE state>
	void operator()(CellPtr& c, const Field* _field) {
		int ix = c->old_data.latx;
		int iy = c->old_data.laty;
		int iz = c->old_data.latz;
		double old_u = c->old_data.ca2p;
		double old_v = c->old_data.ex_inert;
		double old_p = c->old_data.IP3;
		//c->data.ca2p = old_u;
		//c->data.ex_inert = old_v;
		//c->data.IP3 = old_p;
		c->data.ca2p_diff = 0;
		//assert(iz < iz_bound);
		auto& old_a = _field->old_ATP; auto& B = _field->old_ext_stim;
		
		double tmp_A = grid8_avg<double>(old_a, ix, iy, iz);
		double tmp_B = grid8_avg<double>(B, ix, iy, iz);
		c->data.ca2p_diff = fu(c->old_data.ca2p, c->old_data.ex_inert, c->old_data.IP3, tmp_B);
		//new data set same as old
		c->data.ex_inert += DT_Ca*fv<state>(c);
		c->data.IP3 += DT_Ca*fp<state>(tmp_A, c);

		CELL_STATE loop_st[4] = { ALIVE,DEAD,FIX,MUSUME };
		for (int i = 0; i < 4; i++) {
			CELL_STATE& stt = loop_st[i];
			int conn_size = c->old_data.connected_cell[stt].size();
			for (int j = 0; j < conn_size; j++) {
				auto&oppo = c->old_data.connected_cell[stt][j];
				c->data.ca2p_diff  +=
					ca2p_du*IAG<state>(c->old_data.agek)
					*c->old_data.gj[stt][j] * (oppo->old_data.ca2p - c->old_data.ca2p);
				c->data.IP3 += DT_Ca*dp*IAG<state>(c->old_data.agek)
					*c->old_data.gj[stt][j] * (oppo->old_data.IP3 - c->old_data.IP3);

			}
		}
		int conn_size = c->old_data.connected_cell[ALIVE].size();
		for (int i = 0; i < conn_size;i++) {
			auto&oppo = c->old_data.connected_cell[ALIVE][i];
			c->data.gj[ALIVE][i] += DT_Ca*fw(fabs(oppo->old_data.ca2p - c->old_data.ca2p), c->old_data.gj[ALIVE][i]);
		}
		c->data.ca2p += DT_Ca*c->data.ca2p_diff;
		c->data.ca2p_avg += c->data.ca2p;

	}
};

class Field::ca_loop_init {
public:
	template<CELL_STATE state>
	void operator()(CellPtr& c, const Field* _field) {
		Field* field = const_cast<Field*>(_field);
		c->old_data.ca2p = c->data.ca2p;
		c->old_data.ex_inert = c->data.ex_inert;
		c->old_data.IP3 = c->data.IP3;
		c->old_data.ca2p_diff = c->data.ca2p = 0;
		c->old_data.gj= c->data.gj; //copied
	}
};

class Field::ca_result_setter {
public:
	template<CELL_STATE state>
	void operator()(CellPtr& c, const Field* _field) {
		STATIC_IF(state == ALIVE || state == FIX || state == MUSUME) {
			c->data.ca2p_avg /= Ca_ITR;

		}STATIC_ELSE{
			c->data.ca2p_avg = c->data.ca2p = 0;
		}STATIC_ENDIF;
	}
};

void Field::ca_dyn() {
	//air stim is set when map set
	proc_all_state_cell<ca_dyn_init_current,ALIVE,DEAD,FIX,MUSUME,AIR>(a_cell);
	
	int iz_bound = (int)((zzmax + FAC_MAP*R_max) / dz);
	int prev_x, prev_y, prev_z, next_x, next_y, next_z;
	for (int itr = 0; itr < Ca_ITR; itr++) {
		if (itr % 10 == 0)printf("itr:%d\n", itr);
		for (int i = 0; i <= NX; i++) {
			for (int j = 0; j <= NY; j++) {
				for (int k = 0; k <= NZ; k++) {
					old_ATP[i][j][k] = ATP[i][j][k];
				}
			}
		}
		proc_all_state_cell<ca_loop_init, ALIVE, DEAD, FIX, MUSUME, AIR>(a_cell);
		proc_all_state_cell<ca_dead_proc,DEAD>(a_cell);
		proc_all_state_cell<ca_common_cell_proc, ALIVE, FIX, MUSUME>(a_cell);
		for (auto& cdead : a_cell[DEAD]) {
			cdead->old_data.ca2p_diff = cdead->data.ca2p_diff = 0;
		}

		for (int j = 0; j < NX; j++) {
			prev_x= (j - 1 + NX) % NX;
			next_x = (j + 1 + NX) % NX;
			for (int k = 0; k < NY; k++) {
				prev_y = (k - 1 + NY) % NY;
				next_y = (k + 1 + NY) % NY;
				for (int l = 0; l < iz_bound; l++) {
					prev_z = 0, next_z = 0;
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
					auto& old_a = old_ATP; auto& map2 = cell_map2;

					double my_a = old_a[j][k][l];

					CellPtr*& m = cell_map[j][k][l];
					double tmp = m != nullptr ? (*m)->data.ca2p_diff : 0;
					ATP[j][k][l] = old_a[j][k][l]
						+ DT_Ca * (Da * (map2[prev_x][k][l] * (old_a[prev_x][k][l] - my_a)
							+ map2[next_x][k][l] * (old_a[next_x][k][l] - my_a)
							+ map2[j][prev_y][l] * (old_a[j][prev_y][l] - my_a)
							+ map2[j][next_y][l] * (old_a[j][next_y][l] - my_a)
							+ map2[j][k][prev_z] * (old_a[j][k][prev_z] - my_a)
							+ map2[j][k][next_z] * (old_a[j][k][next_z] - my_a)) *inv_dx*inv_dx
							+ fa(tmp, my_a) + air_stim[j][k][l]);
				}
			}

			
		}
		//bc
		for (int l = 0; l < iz_bound; l++) {
			for (int __x=0; __x < NX; __x++)ATP[__x][NY][l] = ATP[__x][0][l];
			for (int __y=0; __y < NY; __y++)ATP[NX][__y][l] = ATP[0][__y][l];
		}
		
	
	}
	
	proc_all_state_cell<ca_result_setter, ALIVE, DEAD, FIX, MUSUME, AIR>(a_cell);
	old_ATP = ATP;
}

void Field::calc_ext_stim() {
	old_ext_stim = ext_stim;
	auto& old_B = old_ext_stim; auto& B = ext_stim;
	int prev_x, prev_y, prev_z, next_x, next_y, next_z;
	int iz_bound = (int)((zzmax + FAC_MAP*R_max) / dz);
	for (int j = 0; j < NX; j++) {
		prev_x = (j - 1 + NX) % NX;
		next_x = (j + 1 + NX) % NX;
		for (int k = 0; k < NY; k++) {
			prev_y = (k - 1 + NY) % NY;
			next_y = (k + 1 + NY) % NY;
			for (int l = 0; l < iz_bound; l++) {
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
				//cond is processed when map setting
				auto& map2 = cell_map2;
				B[j][k][l] = old_B[j][k][l]+
					DT_Ca * (DB * (cell_map2[prev_x][k][l] * (old_B[prev_x][k][l] - old_B[j][k][l])
						+ cell_map2[j][prev_y][l] * (old_B[j][prev_y][l] - old_B[j][k][l])
						+ cell_map2[j][k][prev_z] * (old_B[j][k][prev_z] - old_B[j][k][l])
						+ cell_map2[next_x][k][l] * (-old_B[j][k][l] + old_B[next_x][k][l])
						+ cell_map2[j][next_y][l] * (-old_B[j][k][l] + old_B[j][next_y][l])
						+ cell_map2[j][k][next_z] * (-old_B[j][k][l] + old_B[j][k][next_z])) / (dz * dz)
						+ fB(age_map[j][k][l], old_B[j][k][l], cornif_map[j][k][l]));
			}
		}
	}
	//bc
	for (int l = 0; l <= iz_bound; l++) {
		for (int j = 0; j < NX; j++) B[j][NY][l] = B[j][0][l];
		for (int k = 0; k <= NY; k++) B[NX][k][l] = B[0][k][l];
	}
	old_B = B;
}

