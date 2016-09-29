#include "CellManager.h"

#include <fstream>
#include <cinttypes>
#include <iostream>
CellManager::CellManager(size_t N):Lockfree_push_stack_dyn<Cell*>(N),remove_queue(N)
{
}

Cell* CellManager::create(CELL_STATE _state, int stem_orig_id, double _x, double _y, double _z, double _radius, double _ca2p, double _ca2p_avg, double _IP3, double _ex_inert, double _agek, double _ageb, double _ex_fat, double _in_fat, double _spr_nat_len, int _rest_div_times, bool _is_malignant)
{
    Cell* cptr=new Cell(Cell::ctor_cookie(),
        _state,
        stem_orig_id,
        _ca2p,
        _IP3,
        _ex_inert,
        _agek, _ageb, _ex_fat, _in_fat, _spr_nat_len,
        _x, _y, _z, _rest_div_times,

        _radius, _ca2p_avg, _is_malignant);
    cptr->set_index(this->push_back_with_index(cptr));
    return cptr;
}
void CellManager::add_remove_queue(size_t idx)
{
    remove_queue.push_back(idx);
}

void CellManager::remove_exec()
{
    for (size_t i = 0; i < remove_queue.size(); ++i) {
        size_t remove_idx = remove_queue[i];
        // _data[remove_idx]->removed=true;
        /*
        if (_data[remove_idx]->state == DEAD) {
            _data[remove_idx]->cst.k_disap_timestep = current_timestep;
            cst_store.push_back(std::move(_data[remove_idx]->cst));
        }
        */

        delete _data[remove_idx];
        _data[remove_idx] = _data[--_next];
        _data[remove_idx]->migrate(remove_idx);
    }
    remove_queue.clear();
}

void CellManager::load(const std::string& path) {
    auto&cman = *this;
    /*
    ファイル読み込み試行
    */
    std::ifstream dstrm(path);

    if (!dstrm) {
        throw std::runtime_error("Failed to load the cell data file:"_s + path);
    }

    //using namespace cont;


    /*
    読み込み用一時変数
    */
    CELL_STATE state = UNUSED;
    std::string line;
    int div_times = 0, touch = 0, pair_cell_id = 0, stem_orig_id = 0;
    double rad, ageb, agek, x, y, z, fat, spr_len, ex_fat, ca2p_avg, ca2p;

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
    std::vector<Cell*> tmp_pair(pm->MEMB_NUM_X*pm->MEMB_NUM_Y+100000);

    unsigned int phase = 0;
    int nmemb = 0;
    int nder = 0;

    while (std::getline(dstrm, line)) {

        sscanf(line.c_str(), "%*d %" SCNuFAST8 " %lf %lf %lf %lf %lf %lf %lf %lf %d %lf %lf %d %lf %d %d",
            (uint_fast8_t*)&state, &rad, &ageb, &agek, &ca2p, &x, &y, &z, &ca2p_avg, &div_times, &ex_fat, &fat, &touch, &spr_len, &pair_cell_id, &stem_orig_id);

        /*
        BLANKにたどり着いたら終了
        */
        if (state == BLANK)break;

        /*
        validation
        */
        if (pm->SYSTEM == BASAL && (state == ALIVE || state == DEAD || state == AIR)) {

            throw std::runtime_error(" input date must not contain ALIVE or DEAD in case of BASAL\n");
        }
        if (state == DER && rad != pm->R_der) {
            throw std::runtime_error("radii of DER not consistent with param.h\n");
            
        }
        if (state == MEMB && rad != pm->R_memb) {
            throw std::runtime_error("radii of DER not consistent with param.h\n");
        }
        if (phase == 0 && state != MEMB) {
            if (state != DER) {
                throw std::runtime_error("Wrong cell order:The next block of MEMB must be DER.");
            }
            phase++;
        }

        if (phase == 1 && state != DER) {
            phase++;
        }
        if (phase > 0 && state == MEMB) {
            throw std::runtime_error("Wrong cell order:A MEMB data is found out of MEMB block.");
        }
        if (phase > 1 && state == DER) {
            throw std::runtime_error("Wrong cell order:A DER data is found out of DER block.");
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
            pm->IP3_init,
            pm->ex_inert_init,
            agek, ageb,
            ex_fat, fat,
            spr_len,
            div_times,
            (unsigned int)stem_orig_id < pm->MALIG_NUM
        );

        if(tmp_pair.size()<=id_count)tmp_pair.resize(id_count + 1);
        if (pair_cell_id > -1) {
            if (tmp_pair[id_count] != nullptr) {

                if (tmp_pair[id_count]->pair != nullptr) {
                    throw std::logic_error("Wrong cell pairing");
                }
                cptr->pair = tmp_pair[id_count];
                tmp_pair[id_count]->pair = cptr;
            }
            else {
                tmp_pair.resize(pair_cell_id + 1);
                tmp_pair[pair_cell_id] = cptr;
            }
        }
        std::cout << "Phase " << phase << "  Cell loaded:" << id_count++ << std::endl;

    }
    cman.nmemb = nmemb;
#ifdef TRI_MEMB
    assert(nmemb == (NMX*NMY_tri));
    assert(NMY_tri % 2 == 0);
#endif
    cman.nder = nder;
}