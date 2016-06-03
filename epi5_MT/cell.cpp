#include "cell.h"
#include "cellmanager.h"
#include <iostream>
#include <istream>
#include <fstream>
#include <sstream>
#include <cinttypes>
SwapData<atomic_double[cont::MAX_CELL_NUM][3]> Cell::pos_s;
SwapData<double[cont::MAX_CELL_NUM]> Cell::ca2p_s;
SwapData<double[cont::MAX_CELL_NUM]>Cell::IP3_s;
SwapData<double[cont::MAX_CELL_NUM]>Cell::ex_inert_s;
SwapData<double[cont::MAX_CELL_NUM]>Cell::agek_s;
SwapData<double[cont::MAX_CELL_NUM]>Cell::ageb_s;
SwapData<double[cont::MAX_CELL_NUM]>Cell::ex_fat_s;
SwapData<double[cont::MAX_CELL_NUM]>Cell::in_fat_s;
SwapData<double[cont::MAX_CELL_NUM]>Cell::spr_nat_len_s;
//SwapData<int[cont::MAX_CELL_NUM]>Cell::rest_div_times_s;
SwapData<std::unordered_map<Cell*, double>[cont::MAX_CELL_NUM]>Cell::gj_s;
CellManager Cell::cells;
unsigned int Cell::nmemb = 0;
unsigned int Cell::nder = 0;

void Cell::pos_swap()
{
	pos_s.swap();
}

void Cell::ca2p_swap()
{
	ca2p_s.swap();
}

void Cell::IP3_swap()
{
	IP3_s.swap();
}

void Cell::ex_inert_swap()
{
	ex_inert_s.swap();
}

void Cell::agek_swap()
{
	agek_s.swap();
}

void Cell::ageb_swap()
{
	ageb_s.swap();
}

void Cell::ex_fat_swap()
{
	ex_fat_s.swap();
}

void Cell::in_fat_swap()
{
	in_fat_s.swap();
}

void Cell::spr_nat_len_swap()
{
	spr_nat_len_s.swap();
}

void Cell::load_from_file(std::string path)
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
	CELL_STATE state=UNUSED;
	std::string line;
	int div_times=0, touch=0, pair_cell_id=0, stem_orig_id=0;
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

		auto cptr = Cell::create(
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
	Cell::nmemb = nmemb;
	Cell::nder = nder;
}


size_t Cell::get_index()
{
	return index;
}

void Cell::set_index(size_t i)
{
	index = i;
}

void Cell::migrate(size_t dest_idx)
{
	set_index(dest_idx);
	size_t _index = get_index();
	x._migrate(_index);
	y._migrate(_index);
	z._migrate(_index);
	ca2p._migrate(_index);
	IP3._migrate(_index);
	ex_inert._migrate(_index);
	agek._migrate(_index);
	ageb._migrate(_index);
	in_fat._migrate(_index);
	ex_fat._migrate(_index);
	spr_nat_len._migrate(_index);
	//rest_div_times._migrate(_index);
}

CellPtr Cell::create(CELL_STATE _state, double _x, double _y, double _z, double _radius, double _ca2p, double _ca2p_avg, double _IP3, double _ex_inert, double _agek, double _ageb, double _ex_fat, double _in_fat, double _spr_nat_len, double _div_age_thresh, int _rest_div_times, bool _is_malignant)
{
	CellPtr cptr = std::make_shared<Cell>(_state, _radius, _ca2p_avg, _div_age_thresh, _is_malignant);
	cptr->set_index(Cell::cells.register_cell(cptr));
	size_t _index = cptr->get_index();
	cptr->x.init(_index);
	cptr->y.init(_index);
	cptr->z.init(_index);
	cptr->ca2p.init(_index);
	cptr->IP3.init(_index);
	cptr->ex_inert.init(_index);
	cptr->agek.init(_index);
	cptr->ageb.init(_index);
	cptr->in_fat.init(_index);
	cptr->ex_fat.init(_index);
	cptr->spr_nat_len.init(_index);
	//cptr->rest_div_times.init(_index);

	cptr->x._set(_x);
	cptr->y._set(_y);
	cptr->z._set(_z);
	cptr->ca2p._set(_ca2p);
	cptr->IP3._set(_IP3);
	cptr->ex_inert._set(_ex_inert);
	cptr->agek._set(_agek);
	cptr->ageb._set(_ageb);
	cptr->in_fat._set(_in_fat);
	cptr->ex_fat._set(_ex_fat);
	cptr->spr_nat_len._set(_spr_nat_len);
	cptr->rest_div_times = _rest_div_times;

	return cptr;
}

Cell::Cell(CELL_STATE _state, 
	double _radius , double _ca2p_avg ,
	double _div_age_thresh ,
	bool _is_malignant ) :
	state(_state), radius(_radius), ca2p_avg(_ca2p_avg), div_age_thresh(_div_age_thresh), is_malignant(_is_malignant), diff_u(0) {
}
