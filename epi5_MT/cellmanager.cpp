#include "cellmanager.h"
#include "cell.h"
#include <iostream>
#include <istream>
#include <fstream>
#include <sstream>
#include <cinttypes>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <iomanip>

/**
 *  @file Cellの管理
 */

/**
 *  全ての細胞の位置をアップデート
 */
void pos_copy(CellManager& cman)
{
	cman.all_foreach_parallel_native([](Cell* c) {
		c->x.update();
		c->y.update();
		c->z.update();
	});
}

/**
 *  カルシウム濃度を次の値とスワップ
 */
void ca2p_swap(CellManager& cman)
{
	cman.ca2p_s.swap();
}

/**
 *  IP3濃度を次の値とスワップ
 */
void IP3_swap(CellManager& cman)
{
	cman.IP3_s.swap();
}

/**
 *  角化を行う。
 *  ALIVE->DEAD
 */
void cornificate(CellManager & cman, Cell * const RESTRICT al)
{
	al->state = DEAD;
    al->cst.k_cornified_timestep=cman.current_timestep;
	printf("sw updated:%d\n", ++cman.sw);
}

/**
 *  細胞の登録
 */
size_t CellManager::register_cell(const CellPtr & c)
{
	return push_back_with_index(c);
}

/**
 *  基底膜の初期化
 */
void CellManager::_memb_init()
{
	using namespace cont;
	auto&cells = *this;
	cells.memb_foreach_with_index([&](CellPtr& cptr, size_t j) {
		size_t jj = j%NMX;
		size_t kk = j / NMX;
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
        //assert(cptr->md.memb_l && cptr->md.memb_r && cptr->md.memb_b && cptr->md.memb_u);
		cptr->connected_cell.push_back(cptr->md.memb_l);
		cptr->connected_cell.push_back(cptr->md.memb_r);
		cptr->connected_cell.push_back(cptr->md.memb_b);
		cptr->connected_cell.push_back(cptr->md.memb_u);
	});

    //2pass
    //直接周囲8つ接続することも可能だが、4つ接続するだけで位置関係は定義できているのでそれを利用していきたい

    cells.memb_foreach_with_index([&](CellPtr& cptr, size_t j) {

        cptr->connected_cell.push_back(cptr->md.memb_l->md.memb_u);
        cptr->connected_cell.push_back(cptr->md.memb_b->md.memb_l);
        cptr->connected_cell.push_back(cptr->md.memb_r->md.memb_b);
        cptr->connected_cell.push_back(cptr->md.memb_u->md.memb_r);
    });


    //see cont::MEMB_ADJ_CONN_NUM

}

/**
 *  ファイルから細胞のデータを読み込む
 */
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
    double rad, ageb, agek, x, y, z, fat, spr_len, ex_fat,ca2p_avg,ca2p;

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

        sscanf(line.c_str(), "%*d %" SCNuFAST8 " %lf %lf %lf %lf %lf %lf %lf %lf %d %lf %lf %d %lf %d %d",
            (uint_fast8_t*)&state, &rad, &ageb, &agek,&ca2p, &x, &y, &z, &ca2p_avg,&div_times, &ex_fat, &fat, &touch, &spr_len, &pair_cell_id, &stem_orig_id);

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
            stem_orig_id,
			x, y, z,
			rad,
            ca2p, ca2p_avg,
			IP3_init,
			ex_inert_init,
			agek, ageb,
			ex_fat, fat,
			spr_len,
			div_times,
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

/**
 *  細胞をインデックスで削除キューに追加
 *  (スレッドセーフ)
 */
void CellManager::add_remove_queue(size_t idx)
{
	remove_queue.push_back(idx);
}

/**
 *  削除リストにある細胞を削除
 *  (スレッドセーフではない)
 */
void CellManager::remove_exec()
{
	for (size_t i = 0; i < remove_queue.size(); ++i) {
        size_t remove_idx=remove_queue[i];
       // _data[remove_idx]->removed=true;
        if(_data[remove_idx]->state==DEAD){
        _data[remove_idx]->cst.k_disap_timestep=current_timestep;
        cst_store.push_back(std::move(_data[remove_idx]->cst));
        }
        delete _data[remove_idx];
        _data[remove_idx] = _data[--_next];
        _data[remove_idx]->migrate(remove_idx);
	}
	remove_queue.clear();
}

void CellManager::init_internal(std::string init_data_path)
{
	_load_from_file(init_data_path);
	_memb_init();

}


/**
 *  細胞のインスタンスの生成
 *  (スレッドセーフ)
 */
CellPtr CellManager::create(CELL_STATE _state,int stem_orig_id, double _x, double _y, double _z, double _radius, double _ca2p, double _ca2p_avg, double _IP3, double _ex_inert, double _agek, double _ageb, double _ex_fat, double _in_fat, double _spr_nat_len, int _rest_div_times, bool _is_malignant)
{
    //use smart ptr
    Cell* cptr = new Cell(
                Cell::ctor_cookie(),
                _state,
                stem_orig_id,
                ca2p_s,
                IP3_s,
                _ex_inert,
                _agek,_ageb,_ex_fat,_in_fat,_spr_nat_len,
            _x,_y,_z, _rest_div_times,

                _radius, _ca2p_avg, _is_malignant);

    cptr->set_index(this->register_cell(cptr));
	size_t _index = cptr->get_index();
	cptr->ca2p.init(_index);
	cptr->IP3.init(_index);

	cptr->ca2p._set(_ca2p);
	cptr->IP3._set(_IP3);

    return cptr;
}

/**
 *  周期境界を考慮して全細胞の位置を修正
 */
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

template<typename T>
void bin_write(std::ofstream& ofs,T* ptr){
ofs.write(reinterpret_cast<const char*>(ptr),sizeof(T));
}

template<typename First,typename Second,typename... T>
void bin_write(std::ofstream& ofs,First* fptr,Second* sptr,T*... ptr){
    bin_write<First>(ofs,fptr);
    bin_write<Second,T...>(ofs,sptr,ptr...);
}

/**
 *  局在化の情報を細胞に書き込む
 */
void check_localization(CellManager&cman){
    cman.all_foreach_parallel_native([](Cell*const RESTRICT c){
        bool touch=false;
        if(c->state==ALIVE){
            for(size_t i=0;i<c->connected_cell.size();i++){
                if(c->connected_cell[i]->state==DEAD){
                    touch=true;
                    break;
                }
            }
        }
        c->is_touch=touch;

    });
}

/**
 *  細胞のデータを書き出し
 *  @param [in] filename ファイル名
 *  @param [in] binary_mode バイナリで出力
 */
void CellManager::output(const std::string &filename,bool binary_mode)
{
    std::ofstream wfile;
    if(binary_mode){
    wfile.open(filename,std::ios_base::out|std::ios_base::binary|std::ios_base::trunc);
    }else{
        wfile.open(filename);
    }
    if(!wfile){
        std::cout<<"Output file creation error. Filename:"<<filename<<std::endl;
        exit(1);
    }


    check_localization(*this);

    /*
     *  File Format:
     *
     *  integer         ->signed 32bit
     *  floating point  ->double 64bit (IEEE754)
     *  bool            ->char 8bit
     */


    if(binary_mode){
    this->all_foreach([&](Cell* c){ //no parallel
        //data casting
        double ca2p=(double)(c->ca2p());
        int state=(int)(c->state);
        double x=(double)(c->x());
        double y=(double)(c->y());
        double z=(double)(c->z());
        int index=(int)(c->index);
        double radius=(double)(c->radius);
        double ageb=(double)(c->ageb);
        double agek=(double)(c->agek);
        double ca2p_avg=(double)(c->ca2p_avg);
        int rest_div_times=(int)(c->rest_div_times);
        double ex_fat=(double)(c->ex_fat);
        double in_fat=(double)(c->in_fat);
        char is_touch=c->is_touch?1:0;
        double spr_nat_len=(double)(c->spr_nat_len);
        int pair_index=c->pair==nullptr?(int)-1:(int)(c->pair->index);
        int fix_origin_idx=c->fix_origin;
        bin_write(wfile,
                  &index,
                  &state,
                  &radius,
                  &ageb,
                  &agek,
                  &ca2p,
                  &x,
                  &y,
                  &z,
                  &ca2p_avg,
                  &rest_div_times,
                  &ex_fat,
                  &in_fat,
                  &is_touch,
                  &spr_nat_len,
                  &pair_index,
                  &fix_origin_idx);

    });
    }else{
        using namespace std;
        this->all_foreach([&](Cell* c){
            wfile<<c->index<<" "
                <<c->state<<" "
               <<fixed<<setprecision(15)<<c->radius<<" "
              <<fixed<<setprecision(15)<<c->ageb<<" "
             <<fixed<<setprecision(15)<<c->agek<<" "
                <<fixed<<setprecision(15)<<c->ca2p()<<" "
               <<fixed<<setprecision(15)<<c->x()<<" "
              <<fixed<<setprecision(15)<<c->y()<<" "
             <<fixed<<setprecision(15)<<c->z()<<" "
            <<fixed<<setprecision(15)<<c->ca2p_avg<<" "
               <<c->rest_div_times<<" "
              <<fixed<<setprecision(15)<<c->ex_fat<<" "
             <<fixed<<setprecision(15)<<c->in_fat<<" "
            <<(c->is_touch?1:0)<<" "
               <<fixed<<setprecision(15)<<c->spr_nat_len<<" "
              <<(int)(c->pair==nullptr?(int)-1:(int)(c->pair->index))<<" "
                                                      <<c->fix_origin<<std::endl;
        });
    }
}

/**
 *  色々クリーンアップ
 */
void CellManager::clean_up(){
    fprintf(stdout,"gj cleaning...\n");
    fflush(stdout);
    this->all_foreach_parallel_native([](Cell* const RESTRICT c){
        for(auto it=c->gj.begin();it!=c->gj.end();){
            if(!it->second.is_checked()){
                it=c->gj.erase(it);
            }else{
                ++it;
            }
        }
    });
    fprintf(stdout,"done.\n");
    fflush(stdout);
}

void CellManager::append_stat_data(const std::string &filename){
    if(cst_store.empty())return;
    std::ofstream ofs(filename,std::ios_base::app|std::ios_base::out);
    if(!ofs){std::cout<<"stat data append error"<<std::endl;}
    for(cell_stat& cdata:cst_store){
        ofs<<cdata.k_aging_start_timestep<<" "<<cdata.k_cornified_timestep<<" "<<cdata.k_disap_timestep<<std::endl;
    }
    cst_store.clear();
}
