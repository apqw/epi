#pragma once
#include "cell.h"
#include <tbb/parallel_for.h>
#include <tbb/parallel_for_each.h>
#include <tbb/blocked_range.h>
#include <fstream>
#include <cinttypes>
class CellManager;
#define D_CELL_LOOP_ACCESSOR(prefix,start,end)\
template<class Fn>\
inline void prefix##_foreach_with_index(const Fn& lmbd) {\
	for (size_t i = (start); i < (end); ++i) {\
		lmbd((*this)[i], i);\
	}\
}\
\
template<class Fn>\
inline void prefix##_foreach(const Fn& lmbd) {\
	for (size_t i = (start); i < (end); ++i) {\
		lmbd((*this)[i]);\
	}\
}\
\
template<class Fn>\
inline void prefix##_foreach_parallel(const Fn& lmbd) {\
	tbb::parallel_for<size_t,Fn>((start), (end), lmbd);\
}\
template<class Fn>\
inline void prefix##_foreach_parallel_native(const Fn& lmbd) {\
	\
	tbb::parallel_for_each(&_data[0]+(start),&_data[0]+ (end), lmbd);\
}

struct CellTempStruct {
    CELL_STATE state;
    real rad;
    real ageb, agek, ca2p, x, y, z, ca2p_avg;
    int div_times;
    real ex_fat, fat;
    int touch;
    real spr_len;
    int pair_cell_id, stem_orig_id;
    int id_count;
};

typedef void(*OnCellLoad)(CellManager&, const CellTempStruct&,void*);

struct CellLoadProc {
    static std::vector<Cell*> tmp_pair;
    CellLoadProc();
    void operator()(CellManager&, const CellTempStruct&);
};

class CellManager:public Lockfree_push_stack_dyn<Cell*>
{
    Lockfree_push_stack_dyn<size_t> remove_queue;
    //static std::vector<Cell*> tmp_pair;
    size_t nmemb, nder;
    void memb_init();
    std::atomic<unsigned int> sw=0;
public:
    typedef Lockfree_push_stack_dyn<Cell*> Base;
    friend void cornificate(CellManager& cman, Cell*const RESTRICT c); //in cell_state_renew
    friend CellManager init_gen(int nfix,int der);
    CellManager();
    CellManager(size_t);
    CellManager(const CellManager& cm):Base(cm),nmemb(cm.nmemb),nder(cm.nder),sw(cm.sw.load()),remove_queue(cm.remove_queue){}
    CellManager(CellManager&& cm):Base(std::move(cm)),nmemb(std::move(cm.nmemb)),nder(std::move(cm.nder)),sw(cm.sw.load()),remove_queue(std::move(cm.remove_queue)){}
    CellManager& operator=(const CellManager& cm){
    	Base::operator=(cm);
    	nmemb=cm.nmemb;
    	nder=cm.nder;
    	sw=cm.sw.load();
    	remove_queue=cm.remove_queue;
        return *this;
    }

    CellManager& operator=(CellManager&& cm){
        	Base::operator=(cm);
        	nmemb=std::move(cm.nmemb);
        	nder=std::move(cm.nder);
        	sw=cm.sw.load();
        	remove_queue=std::move(cm.remove_queue);
            return *this;
        }
    Cell* create(CELL_STATE _state, int stem_orig_id, real _x = 0, real _y = 0, real _z = 0,
        real _radius = pm->R_max, real _ca2p = pm->ca2p_init, real _ca2p_avg = pm->ca2p_init,
        real _IP3 = pm->IP3_init, real _ex_inert = pm->ex_inert_init,
        real _agek = 0, real _ageb = 0,
        real _ex_fat = 0, real _in_fat = 0,
        real _spr_nat_len = 0,
        int _rest_div_times = 0,
        bool _is_malignant = false);

    Cell* create_resizable(CELL_STATE _state, int stem_orig_id, real _x = 0, real _y = 0, real _z = 0,
            real _radius = pm->R_max, real _ca2p = pm->ca2p_init, real _ca2p_avg = pm->ca2p_init,
            real _IP3 = pm->IP3_init, real _ex_inert = pm->ex_inert_init,
            real _agek = 0, real _ageb = 0,
            real _ex_fat = 0, real _in_fat = 0,
            real _spr_nat_len = 0,
            int _rest_div_times = 0,
            bool _is_malignant = false);

    void add_remove_queue(size_t idx);
    void remove_exec();
    void output(const std::string&filename,bool binary_mode=false);
    void init_value();

    template<typename Fn=CellLoadProc>
    void load(const std::string& path, Fn on=CellLoadProc()) {
        auto&cman = *this;
        /*
        �t�@�C���ǂݍ��ݎ��s
        */
        std::ifstream dstrm(path);

        if (!dstrm) {
            throw std::runtime_error("Failed to load the cell data file:"_s + path);
        }


        /*
        �ǂݍ��ݗp�ꎞ�ϐ�
        */
        CELL_STATE state = UNUSED;
        std::string line;
        int div_times = 0, touch = 0, pair_cell_id = 0, stem_orig_id = 0;
        double rad, ageb, agek, x, y, z, fat, spr_len, ex_fat, ca2p_avg, ca2p;

        /*
        �זE�̃C���f�b�N�X�̃J�E���g
        (1�s��1�זE����̂ŁA1�s�ǂݍ��ނ��ƂɃC���N�������g����)
        */
        int id_count = 0;

        /*
        �y�A���ꎞ�I�ɕۑ�����vector

        �y�A�����זE�����߂Č������Ƃ��A���̃y�A�̃C���f�b�N�X�͕K����������Ȃ̂ł܂���������Ă��Ȃ��B
        ���̂��߁A�y�A�����זE������΁A�����̃C���f�b�N�X���L�[�Ƃ��Ĉꎞ�I�ɓ��Ă����A���̃y�A�ɂ��ǂ蒅�����Ƃ��Ɍ݂��̃I�u�W�F�N�g�ɓo�^����B
        */
        //tmp_pair.resize(pm->MEMB_NUM_X*pm->MEMB_NUM_Y + 100000, nullptr);
        //std::fill(tmp_pair.begin(), tmp_pair.end(), nullptr);

        unsigned int phase = 0;
        int nmemb = 0;
        int nder = 0;

        while (std::getline(dstrm, line)) {

            sscanf(line.c_str(), "%*d %" SCNuFAST8 " %lf %lf %lf %lf %lf %lf %lf %lf %d %lf %lf %d %lf %d %d",
                (uint_fast8_t*)&state, &rad, &ageb, &agek, &ca2p, &x, &y, &z, &ca2p_avg, &div_times, &ex_fat, &fat, &touch, &spr_len, &pair_cell_id, &stem_orig_id);

            /*
            BLANK�ɂ��ǂ蒅������I��
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
            CellTempStruct cts;
            cts.ageb = (real)ageb; cts.agek = (real)agek; cts.ca2p = (real)ca2p; cts.ca2p_avg = (real)ca2p_avg;
            cts.div_times = div_times; cts.ex_fat = (real)ex_fat; cts.fat = (real)fat;
            cts.pair_cell_id = pair_cell_id; cts.rad = (real)rad; cts.spr_len = (real)spr_len;
            cts.state = state; cts.stem_orig_id = stem_orig_id; cts.touch = touch;
            cts.x = (real)x; cts.y = (real)y; cts.z = (real)z; cts.id_count = id_count;
            on(cman, cts);

            //std::cout << "Phase " << phase << "  Cell loaded:" << id_count++ << std::endl;
            id_count++;

        }
        cman.nmemb = nmemb;
        if(nmemb!=pm->MEMB_NUM_X*pm->MEMB_NUM_Y){
        	throw std::runtime_error("MEMB num of input file is not same with this.");
        }
        cman.nder = nder;
    }

    D_CELL_LOOP_ACCESSOR(all, 0, size());
    D_CELL_LOOP_ACCESSOR(memb, 0, nmemb);
    D_CELL_LOOP_ACCESSOR(der, nmemb, nder + nmemb);
    D_CELL_LOOP_ACCESSOR(other, nder + nmemb, size());
    D_CELL_LOOP_ACCESSOR(non_memb, nmemb, size());
    void _memb_init();
    void pos_update();
    void pos_periodic_fix();
    bool should_calc_ca();
    void ca_calc_condition_reset();
    //~CellManager();
};

