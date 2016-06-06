#include "cellmanager.h"
#include "cell.h"
#include <iostream>
#include <istream>
#include <fstream>
#include <sstream>
#include <cinttypes>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
void CellManager::pos_swap()
{
	//pos_s.swap();
}

void CellManager::pos_copy()
{
	tbb::parallel_for(tbb::blocked_range<size_t>(0, size()), [&](const tbb::blocked_range<size_t>& range) {
		for (size_t i = range.begin(); i != range.end(); ++i) {
			(*this)[i]->x.update();
			(*this)[i]->y.update();
			(*this)[i]->z.update();
		}
	});
}

void CellManager::ca2p_swap()
{
	ca2p_s.swap();
}

void CellManager::IP3_swap()
{
	IP3_s.swap();
}

void CellManager::ex_inert_swap()
{
	ex_inert_s.swap();
}

void CellManager::agek_swap()
{
	//agek_s.swap();
}

void CellManager::ageb_swap()
{
	//ageb_s.swap();
}

void CellManager::ex_fat_swap()
{
	//ex_fat_s.swap();
}

void CellManager::in_fat_swap()
{
	//in_fat_s.swap();
}

void CellManager::spr_nat_len_swap()
{
	//spr_nat_len_s.swap();
}

size_t CellManager::register_cell(const CellPtr & c)
{
	return push_back_with_index(c);
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

CellPtr CellManager::create(CELL_STATE _state, double _x, double _y, double _z, double _radius, double _ca2p, double _ca2p_avg, double _IP3, double _ex_inert, double _agek, double _ageb, double _ex_fat, double _in_fat, double _spr_nat_len, double _div_age_thresh, int _rest_div_times, bool _is_malignant)
{
	
	std::shared_ptr<Cell> cptr = std::make_shared<Cell>(
		Cell::ctor_cookie(),
		_state,
		ca2p_s,
		IP3_s,
		ex_inert_s,
		_agek,_ageb,_ex_fat,_in_fat,_spr_nat_len,
	_x,_y,_z,
		
		_radius, _ca2p_avg, _div_age_thresh, _is_malignant);
	cell_store.push_back(cptr);
	cptr->set_index(this->register_cell(cptr.get()));
	size_t _index = cptr->get_index();
	cptr->ca2p.init(_index);
	cptr->IP3.init(_index);
	cptr->ex_inert.init(_index);
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
	cptr->ex_inert._set(_ex_inert);
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

void cman_load_from_file(CellManager & cman, std::string path)
{
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
			&state, &rad, &ageb, &agek, &x, &y, &z, &div_times, &ex_fat, &fat, &touch, &spr_len, &pair_cell_id, &stem_orig_id);

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

void cell_pos_periodic_fix(CellManager& cman) {
	
	cman.all_foreach_parallel_native([&](Cell* c) {
		using namespace cont;
		if (c->x._next() > LX) {
			c->x._next() -= LX;
		}
		else if (c->x._next() < 0) {
			c->x._next() += LX;
		}

		if (c->y._next() > LY) {
			c->y._next() -= LY;
		}
		else if (c->y._next() < 0) {
			c->y._next() += LY;
		}
	});
	
}