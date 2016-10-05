#include "CellManager.h"

#include <fstream>
#include <cinttypes>
#include <iostream>
#include <iomanip>
std::vector<Cell*> CellLoadProc::tmp_pair;
CellManager::CellManager(size_t N):Lockfree_push_stack_dyn<Cell*>(N),remove_queue(N)
{
}
CellManager::CellManager() : Lockfree_push_stack_dyn<Cell*>(0), remove_queue(0)
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
//spelling: resizeable?
Cell* CellManager::create_resizable(CELL_STATE _state, int stem_orig_id, double _x, double _y, double _z, double _radius, double _ca2p, double _ca2p_avg, double _IP3, double _ex_inert, double _agek, double _ageb, double _ex_fat, double _in_fat, double _spr_nat_len, int _rest_div_times, bool _is_malignant)
{
	test_realloc();
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

        delete _data[remove_idx];
        _data[remove_idx] = _data[--_next];
        _data[remove_idx]->migrate(remove_idx);
    }
    remove_queue.clear();
}

CellLoadProc::CellLoadProc() {
    tmp_pair.resize(pm->MEMB_NUM_X*pm->MEMB_NUM_Y + 100000, nullptr);
}

void CellLoadProc::operator()(CellManager&cman, const CellTempStruct& cts) {
    //cman.test_realloc();
    auto cptr = cman.create_resizable(
        cts.state,
        cts.stem_orig_id,
        cts.x, cts.y, cts.z,
        cts.rad,
        cts.ca2p, cts.ca2p_avg,
        pm->IP3_init,
        pm->ex_inert_init,
        cts.agek, cts.ageb,
        cts.ex_fat, cts.fat,
        cts.spr_len,
        cts.div_times,
        (unsigned int)cts.stem_orig_id < pm->MALIG_NUM
    );

    if (tmp_pair.size() <= cts.id_count)tmp_pair.resize(cts.id_count + 1);
    if (cts.pair_cell_id > -1) {
        if (tmp_pair[cts.id_count] != nullptr) {

            if (tmp_pair[cts.id_count]->pair != nullptr) {
                throw std::logic_error("Wrong cell pairing");
            }
            cptr->pair = tmp_pair[cts.id_count];
            tmp_pair[cts.id_count]->pair = cptr;
        }
        else {
            tmp_pair.resize(cts.pair_cell_id + 1);
            tmp_pair[cts.pair_cell_id] = cptr;
        }
    }
    
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
        double ca2p=(double)(c->ca2p);
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
                <<fixed<<setprecision(15)<<c->ca2p<<" "
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

void CellManager::init_value(){
	all_foreach_parallel_native([](Cell*c){
		c->ca2p=pm->ca2p_init;
		c->ca2p_avg=pm->ca2p_init;
		c->ex_inert=pm->ex_inert_init;
		c->IP3=pm->IP3_init;
	});
	other_foreach_parallel_native([](Cell*c){
		if(c->state==DEAD){
			c->ca2p_avg=0.0;
			c->ca2p=0.0;
		}
	});
}


/**
*  基底膜の初期化
*/
void CellManager::_memb_init()
{
    auto&cells = *this;
    if (pm->USE_TRI_MEMB) {
        //1pass
        cells.memb_foreach_with_index([&](Cell* cptr, size_t j) {
            cptr->md.memb.resize(6);
            size_t jj = j%pm->MEMB_NUM_X;
            size_t kk = j / pm->MEMB_NUM_X;
            if (jj == 0) { //most left
                cptr->md.memb[3] = cells[j + pm->MEMB_NUM_X - 1];
            }
            else {
                cptr->md.memb[3] = cells[j - 1];
            }

            if (jj == pm->MEMB_NUM_X - 1) {
                cptr->md.memb[0] = cells[j - (pm->MEMB_NUM_X - 1)];
            }
            else {
                cptr->md.memb[0] = cells[j + 1];
            }

        });
        //2pass
        cells.memb_foreach_with_index([&](Cell* cptr, size_t j) {

            size_t jj = j%pm->MEMB_NUM_X;
            size_t kk = j / pm->MEMB_NUM_X;
            // size_t YNUM=NMY/y_comp_ratio;
            size_t top = (j - pm->MEMB_NUM_X + cells.nmemb) % (cells.nmemb);
            size_t bot = (j + pm->MEMB_NUM_X) % (cells.nmemb);
            if (kk % 2 == 0) {

                cptr->md.memb[1] = cells[top];
                cptr->md.memb[2] = cells[top]->md.memb[3];
                cptr->md.memb[4] = cells[bot]->md.memb[3];
                cptr->md.memb[5] = cells[bot];

            }
            else {
                cptr->md.memb[1] = cells[top]->md.memb[0];
                cptr->md.memb[2] = cells[top];
                cptr->md.memb[4] = cells[bot];
                cptr->md.memb[5] = cells[bot]->md.memb[0];
            }
            for (int i = 0; i < 6; i++) {
                cptr->connected_cell.push_back(cptr->md.memb[i]);
            }
        });
    }
    else {
        cells.memb_foreach_with_index([&](Cell* cptr, size_t j) {
            cptr->md.memb.resize(4);
            size_t jj = j%pm->MEMB_NUM_X;
            size_t kk = j / pm->MEMB_NUM_X;
            if (jj == 0) {
                cptr->md.memb[2] = cells[j + pm->MEMB_NUM_X - 1];
            }
            else {
                cptr->md.memb[2] = cells[j - 1];
            }


            if (jj == pm->MEMB_NUM_X - 1) {
                cptr->md.memb[0]= cells[j - (pm->MEMB_NUM_X - 1)];
            }
            else {
                cptr->md.memb[0] = cells[j + 1];
            }

            if (kk == 0) {
                cptr->md.memb[3] = cells[j + pm->MEMB_NUM_X*pm->MEMB_NUM_Y - pm->MEMB_NUM_X];
            }
            else {
                cptr->md.memb[3] = cells[j - pm->MEMB_NUM_X];
            }


            if (kk == pm->MEMB_NUM_Y - 1) {
                cptr->md.memb[1] = cells[j - (pm->MEMB_NUM_X*pm->MEMB_NUM_Y - pm->MEMB_NUM_X)];
            }
            else {
                cptr->md.memb[1] = cells[j + pm->MEMB_NUM_X];
            }
            //assert(cptr->md.memb_l && cptr->md.memb_r && cptr->md.memb_b && cptr->md.memb_u);
            for (int i = 0; i < 4; i++) {
                cptr->connected_cell.push_back(cptr->md.memb[i]);
            }

        });
    }


    //see cont::MEMB_ADJ_CONN_NUM
    //cell_shuffle(0, nmemb - 1);
}

void CellManager::pos_update() {
    this->all_foreach_parallel_native([](Cell* c) {
        c->x.update();
        c->y.update();
        c->z.update();
    });
}