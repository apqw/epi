#include "cellmanager.h"
#include "cell.h"
#include <iostream>
#include <istream>
#include <fstream>
#include <sstream>
#include <cinttypes>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
void pos_copy(CellManager& cman)
{
	cman.all_foreach_parallel_native([](Cell* c) {
		c->x.update();
		c->y.update();
		c->z.update();
	});
}

void ca2p_swap(CellManager& cman)
{
	cman.ca2p_s.swap();
}

void IP3_swap(CellManager& cman)
{
	cman.IP3_s.swap();
}

void cornificate(CellManager & cman, Cell * const RESTRICT al)
{
	al->state = DEAD;
	printf("sw updated:%d\n", ++cman.sw);
}


size_t CellManager::register_cell(const CellPtr & c)
{
	return push_back_with_index(c);
}


void CellManager::_memb_init()
{
	using namespace cont;
	auto&cells = *this;
	cells.memb_foreach_with_index([&](CellPtr& cptr, size_t j) {
		int jj = j%NMX;
		int kk = j / NMX;
		if (jj == 0) {
			cptr->md.memb_l = cells[j + NMX - 1];
		}
		else {
			cptr->md.memb_l = cells[j - 1];
		}

		if (jj <= 1) {
			cptr->md.memb_ll = cells[j + NMX - 2];
		}
		else {
			cptr->md.memb_ll = cells[j - 2];
		}

		if (jj == NMX - 1) {
			cptr->md.memb_r = cells[j - (NMX - 1)];
		}
		else {
			cptr->md.memb_r = cells[j + 1];
		}

		if (kk == 0) {
			cptr->md.memb_b = cells[j + NMX*NMY - NMX];
		}
		else {
			cptr->md.memb_b = cells[j - NMX];
		}

		if (kk <= 1) {
			cptr->md.memb_bb = cells[j + NMX*NMY - 2 * NMX];
		}
		else {
			cptr->md.memb_bb = cells[j - 2 * NMX];
		}

		if (kk == NMY - 1) {
			cptr->md.memb_u = cells[j - (NMX*NMY - NMX)];
		}
		else {
			cptr->md.memb_u = cells[j + NMX];
		}

		cptr->connected_cell.push_back(cptr->md.memb_l);
		cptr->connected_cell.push_back(cptr->md.memb_r);
		cptr->connected_cell.push_back(cptr->md.memb_b);
		cptr->connected_cell.push_back(cptr->md.memb_u);
		//test
		/*
		cptr->gj._emplace(cptr->md.memb_l, 0);
		cptr->gj._emplace(cptr->md.memb_r, 0);
		cptr->gj._emplace(cptr->md.memb_b, 0);
		cptr->gj._emplace(cptr->md.memb_u, 0);
		*/
	});
}

void CellManager::_load_from_file(std::string path)
{
	auto&cman = *this;
	/*
	ファイル読み込み試行
	*/
	std::ifstream dstrm;
	dstrm.open(path, std::ios::in);

	if (!dstrm) {
		assert(dstrm);
		throw "load failed";
	}

	using namespace cont;


	/*
	読み込み用一時変数
	*/
	CELL_STATE state = UNUSED;
	std::string line;
	int div_times = 0, touch = 0, pair_cell_id = 0, stem_orig_id = 0;
	double rad, ageb, agek, x, y, z, fat, spr_len, ex_fat;

	/*
	細胞のインデックスのカウント
	(1行に1細胞あるので、1行読み込むごとにインクリメントする)
	*/
	int id_count = 0;

	/*
	ペアを一時的に保存するvector

	ペアを持つ細胞を初めて見つけたとき、そのペアのインデックスは必ず自分より後なのでまだ生成されていない。
	そのため、ペアを持つ細胞があれば、自分のインデックスをキーとして一時的に入れておき、そのペアにたどり着いたときに互いのオブジェクトに登録する。
	*/
	std::vector<CellPtr> tmp_pair(MAX_CELL_NUM, nullptr);

	unsigned int phase = 0;
	int nmemb = 0;
	int nder = 0;

	while (std::getline(dstrm, line)) {

		sscanf(line.c_str(), "%*d %" SCNuFAST8 " %lf %lf %lf %*f %lf %lf %lf %*f %d %lf %lf %d %lf %d %d",
            (uint_fast8_t*)&state, &rad, &ageb, &agek, &x, &y, &z, &div_times, &ex_fat, &fat, &touch, &spr_len, &pair_cell_id, &stem_orig_id);

		/*
		BLANKにたどり着いたら終了
		*/
		if (state == BLANK)break;

		/*
		validation
		*/
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

		auto cptr = cman.create(
			state,
			x, y, z,
			rad,
			ca2p_init, ca2p_init,
			IP3_init,
			ex_inert_init,
			agek, ageb,
			ex_fat, fat,
			spr_len,
			state == FIX ? agki_max_fix : state == MUSUME ? agki_max : 0,
			state == FIX ? div_max : div_times,
			stem_orig_id < MALIG_NUM
			);

		if (pair_cell_id > -1) {
			if (tmp_pair[id_count] != nullptr) {
				assert(tmp_pair[id_count]->pair == nullptr);
				cptr->pair = tmp_pair[id_count];
				tmp_pair[id_count]->pair = cptr;
			}
			else {
				tmp_pair[pair_cell_id] = cptr;
			}
		}
		printf("Phase %d  Cell loaded:%d\n", phase, id_count++);

	}
	cman.nmemb = nmemb;
	cman.nder = nder;
}

void CellManager::add_remove_queue(size_t idx)
{
	remove_queue.push_back(idx);
}

void CellManager::remove_exec()
{
	for (size_t i = 0; i < remove_queue.size(); ++i) {
		_data[i] = _data[--_next];
		_data[i]->migrate(i);
	}
	remove_queue.clear();
}

void CellManager::init_internal(std::string init_data_path)
{
	_load_from_file(init_data_path);
	_memb_init();

}

CellPtr CellManager::create(CELL_STATE _state, double _x, double _y, double _z, double _radius, double _ca2p, double _ca2p_avg, double _IP3, double _ex_inert, double _agek, double _ageb, double _ex_fat, double _in_fat, double _spr_nat_len, double _div_age_thresh, int _rest_div_times, bool _is_malignant)
{
	
	std::shared_ptr<Cell> cptr = std::make_shared<Cell>(
		Cell::ctor_cookie(),
		_state,
		ca2p_s,
		IP3_s,
		_ex_inert,
		_agek,_ageb,_ex_fat,_in_fat,_spr_nat_len,
	_x,_y,_z,
		
		_radius, _ca2p_avg, _div_age_thresh, _is_malignant);
	cell_store.push_back(cptr);
	cptr->set_index(this->register_cell(cptr.get()));
	size_t _index = cptr->get_index();
	cptr->ca2p.init(_index);
	cptr->IP3.init(_index);
	//cptr->ex_inert.init(_index);
	/*
	cptr->agek.init(_index);
	cptr->ageb.init(_index);
	cptr->in_fat.init(_index);
	cptr->ex_fat.init(_index);
	cptr->spr_nat_len.init(_index);
	*/
	//cptr->rest_div_times.init(_index);

	cptr->ca2p._set(_ca2p);
	cptr->IP3._set(_IP3);
	//cptr->ex_inert._set(_ex_inert);
	/*
	cptr->agek._set(_agek);
	cptr->ageb._set(_ageb);
	cptr->in_fat._set(_in_fat);
	cptr->ex_fat._set(_ex_fat);
	cptr->spr_nat_len._set(_spr_nat_len);
	*/
	cptr->rest_div_times = _rest_div_times;

	return cptr.get();
}


void cell_pos_periodic_fix(CellManager& cman) {
	
	cman.all_foreach_parallel_native([&](Cell*const c) {
		using namespace cont;
		if (c->x._next() > LX) {
			c->x._set(c->x._next()-LX);
		}
		else if (c->x._next() < 0) {
			c->x._set(c->x._next() + LX);
		}

		if (c->y._next() > LY) {
			c->y._set(c->y._next() - LY);
		}
		else if (c->y._next() < 0) {
			c->y._set(c->y._next() + LY);
		}
	});
	
}

