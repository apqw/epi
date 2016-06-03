#ifndef CELL_H
#define CELL_H
#include "atomics.h"
#include "define.h"
#include <memory>
#include <vector>
#include <unordered_map>
using CellPtr = std::shared_ptr<Cell>;



class Cell
{
    template<T& sw>
    struct SwapArrAccessor{
        size_t idx;
        const auto& operator()(){
            return sw.first()[idx];
        }
        template<typename U>
        SwapArrAccessor& operator+=(const U& v){
            sw.second()[idx]+=v;
        }
template<typename U>
        void _set(U& v){
sw.first()[idx]=v;
        }

        SwapArrAccessor(const size_t _idx):idx(_idx){}
        void _migrate(const size_t dest_idx){
            sw.first()[dest_idx]=sw.first()[idx];
            sw.second()[dest_idx]=sw.second()[idx];
            idx=dest_idx;
        }
    };

    template<T& sw,size_t subidx>
    struct SwapArrAccessor{
        const size_t idx;
        const auto& operator()(){
            return sw.first()[idx][subidx];
        }
        template<typename U>
        SwapArrAccessor& operator+=(const U& v){
            sw.second()[idx][subidx]+=v;
        }

        template<typename U>
                void _set(U& v){
        sw.first()[idx][subidx]=v;
                }

        SwapArrAccessor(const size_t _idx):idx(_idx){}
        void _migrate(const size_t dest_idx){
            sw.first()[dest_idx][subidx]=sw.first()[idx][subidx];
            sw.second()[dest_idx][subidx]=sw.second()[idx][subidx];
            idx=dest_idx;
        }
    };

    static SwapData<atomic_double[cont::MAX_CELL_NUM][3]> pos_s;
    static SwapData<double[cont::MAX_CELL_NUM]>ca2p_s;
    static SwapData<double[cont::MAX_CELL_NUM]>IP3_s;
    static SwapData<double[cont::MAX_CELL_NUM]>ex_inert_s;
    static SwapData<double[cont::MAX_CELL_NUM]>agek_s;
    static SwapData<double[cont::MAX_CELL_NUM]>ageb_s;
    static SwapData<double[cont::MAX_CELL_NUM]>exfat_s;
    static SwapData<double[cont::MAX_CELL_NUM]>infat_s;
    static SwapData<double[cont::MAX_CELL_NUM]>spr_nat_len_s;
    static SwapData<int[cont::MAX_CELL_NUM]>rest_div_times_s;
    static SwapData<std::unordered_map<Cell*,double>[cont::MAX_CELL_NUM]>gj_s;
    static
public:
    /*
     * 位置
     * 並列で加算など可能
     * 相互作用の計算で同時に同じ細胞の値に加算する可能性がある
     */
    SwapArrAccessor<pos_s,0> x;
    SwapArrAccessor<pos_s,1> y;
    SwapArrAccessor<pos_s,2> z;


    /*
     * 接続された細胞（へのポインタ）のリスト
     * 追加のみ並列実行できる
     * タイムステップを越えてポインタを保持しなく、大量に追加することが考えられるので、コストの安い生のポインタを使用
     */

    Lockfree_push_stack<Cell*,cont::MAX_CONNECTED_CELL_NUM> connected_cell;

    /*
     * 分裂中のペアへのポインタ
     */
    CellPtr pair;

    /*
     * 半径
     * 書き換えないので並列への対応は無し
     */
    const double radius;

    /*
     * Ca2+濃度(平均化されていない)
     * この値はCa2+の計算部分以外では使われないはず
     */
    SwapArrAccessor<ca2p_s> ca2p;

    /*
     * Ca2+濃度(平均)
     */
    double ca2p_avg;

    /*
     * IP3
     */
    SwapArrAccessor<IP3_s> IP3;

    /*
     * 不活性化物質
     */
    SwapArrAccessor<ex_inert_s> ex_inert;

    /*
     * 分化後の年齢
     */
    SwapArrAccessor<agek_s> agek;

    /*
     * 分化前の年齢
     */
    SwapArrAccessor<ageb_s> ageb;

    /*
     * 細胞外脂質
     */
    SwapArrAccessor<ex_fat_s> ex_fat;

    /*
     * 細胞内脂質
     */
    SwapArrAccessor<in_fat_s> in_fat;

    /*
     * 細胞分裂中のバネの自然長
     *
     * 他の細胞の値も書き換えるが、同時に書き換えないようにできる。
     */

    SwapArrAccessor<spr_nat_len_s> spr_nat_len;

    /*
     * 分裂開始年齢のしきい値
     *
     * この値に確率に関する係数がかかって分裂開始年齢になる
     */
    const double div_age_thresh;

    /*
     * 残りの分裂回数
     *
     * 同時に書き換えないようにできる
     */
    SwapArrAccessor<rest_div_times_s> rest_div_times;

    /*
     * 悪性かどうか
     */
    bool is_malignant;

    /*
     * 格子としての位置
     */
    int lat[3];

    /*
     * gj
     */
    SwapArrAccessor<gj_s> gj;

    /*
     * diff_u
     * Ca2+のみ
     */
    double diff_u;

    const size_t index;

    CELL_STATE state;

    Cell(CELL_STATE _state,double _x=0,double _y=0,double _z=0,
         double _radius=cont::R_max,double _ca2p=cont::ca2p_init,double _ca2p_avg = cont::ca2p_init,
         double _IP3=cont::IP3_init,double _ex_inert=cont::ex_inert_init,
         double _agek=0,double _ageb=0,
         double _ex_fat=0,double _in_fat=0,
         double _spr_nat_len=0,
         double _div_age_thresh=0,
         double _rest_div_times=0,
         bool _is_malignant=false,size_t _index=Cell::unused_id[Cell::id_head++]):
    state(_state),index(_index),radius(_radius),ca2p_avg(_ca2p_avg),div_age_thresh(_div_age_thresh),x(_index),y(_index),z(_index),ca2p(_index),IP3(_index),ex_inert(_index),agek(_index),ageb(_index)
    ,ex_fat(_index),in_fat(_index),spr_nat_len(_index),rest_div_times(_index),gj(_index),is_malignant(_is_malignant),diff_u(0){
x._set(_x);
y._set(_y);
z._set(_z);
ca2p._set(_ca2p);
IP3._set(_IP3);
ex_inert._set(_ex_inert);
agek._set(_agek);
ageb._set(_ageb);
in_fat._set(_in_fat);
ex_fat._set(_ex_fat);
spr_nat_len._set(_spr_nat_len);
rest_div_times._set(_rest_div_times);
    }


};

class CellMan {
    std::vector<CellPtr> cell;
    Lockfree_push_stack<int,cont::MAX_CELL_NUM> remove_index;
    //std::vector<int> remove_index;
    Lockfree_push_stack<CellPtr,cont::MAX_CELL_NUM> add_cell;
    //std::vector<CellPtr> add_cell;
    int memb_num; int der_num;
    bool nmemb_is_set; int nder_is_set;
public:
    void set_memb_num(int n) {
        assert(!nmemb_is_set);
        memb_num = n;
        nmemb_is_set = true;
    }

    void set_der_num(int n) {
        assert(!nder_is_set);
        der_num = n;
        nder_is_set = true;
    }


    CellMan():nmemb_is_set(false),nder_is_set(false) {
        cell.reserve(cont::MAX_CELL_NUM);
    }

    void add_direct(CellPtr& made_ptr) {
        cell.push_back(made_ptr);
    }
    void add_queue(const CellPtr& made_ptr) {
        add_cell.push_back(made_ptr);
    }

    void remove_queue(int idx) {
        remove_index.push_back(idx);
    }

    std::vector<CellPtr>& _raw_cell_set() {
        return cell;
    }



    /*
    remove instantly
    */
    template<class L>
    void remove_if(const L& pred) {
        auto it = std::remove_if(cell.begin(), cell.end(), pred);
        cell.erase(it, cell.end());
    }

    template<class L>
    void remove_if_queue(const L& pred) {
        size_t sz = cell.size();
        for (int i = 0; i < sz; ++i) {
            if (pred(cell[i], i)) {
                remove_queue(i);
            }
        }
    }

    template<class L>
    void foreach(const L& lmbd) {
        size_t sz = cell.size();

        for (int i = 0; i < sz; ++i) {
            lmbd(cell[i], i);
        }

    }

    template<class L>
    void foreach_parallel(const L& lmbd) {
        size_t sz = cell.size();
        tbb::parallel_for(tbb::blocked_range<size_t>(0, sz), [&](const tbb::blocked_range< size_t >& range) {
            for (size_t i = range.begin(); i != range.end(); ++i) {
                lmbd(cell[i], i);
            }
        });
    }

    template<class L>
    void foreach_parallel_native(const L& lmbd) {
        tbb::parallel_for_each(cell.begin(), cell.end(), lmbd);
    }




    template<class L>
    void memb_foreach(const L& lmbd) {
        assert(nmemb_is_set);


        for (int i = 0; i < memb_num; ++i) {
        lmbd(cell[i], i);
        }

    }

    template<class L>
    void memb_foreach_parallel(const L& lmbd) {
        assert(nmemb_is_set);

        tbb::parallel_for(tbb::blocked_range<int>(0, memb_num), [&](const tbb::blocked_range< int >& range) {
            for (int i = range.begin(); i != range.end(); ++i) {
                lmbd(cell[i], i);
            }
        });
    }

    template<class L>
    void memb_foreach_parallel_native(const L& lmbd) {
        assert(nmemb_is_set);
        tbb::parallel_for_each(cell.begin(), cell.begin() + memb_num, lmbd);
        /*
        tbb::parallel_for(tbb::blocked_range<int>(0, memb_num), [&](const tbb::blocked_range< int >& range) {
            for (int i = range.begin(); i != range.end(); ++i) {
                lmbd(cell[i], i);
            }
        });
        */
    }

    template<class L>
    void non_memb_foreach_parallel_native(const L& lmbd) {
        assert(nmemb_is_set);
        tbb::parallel_for_each(cell.begin()+memb_num, cell.end(), lmbd);
        /*
        tbb::parallel_for(tbb::blocked_range<int>(0, memb_num), [&](const tbb::blocked_range< int >& range) {
        for (int i = range.begin(); i != range.end(); ++i) {
        lmbd(cell[i], i);
        }
        });
        */
    }

    template<class L>
    void der_foreach(const L& lmbd) {
        assert(nder_is_set);
        for (int i = memb_num + 1; i < memb_num + der_num; ++i) {
            lmbd(cell[i], i);
        }
    }

    template<class L>
    void other_foreach(const L& lmbd) {
        assert(nmemb_is_set);
        assert(nder_is_set);
        size_t sz = cell.size();
        for (int i = memb_num +der_num; i < sz; ++i) {
            lmbd(cell[i], i);
        }
    }

    template<class L>
    void other_foreach_parallel(const L& lmbd) {
        assert(nmemb_is_set);
        assert(nder_is_set);
        size_t sz = cell.size();
        tbb::parallel_for(tbb::blocked_range<int>(memb_num + der_num , sz), [&](const tbb::blocked_range< int >& range) {
            for (int i = range.begin(); i != range.end(); ++i) {
                lmbd(cell[i], i);
            }
        });
    }

    template<class L>
    void other_foreach_parallel_native(const L& lmbd) {
        assert(nmemb_is_set);
        assert(nder_is_set);
        tbb::parallel_for_each(cell.begin() + memb_num + der_num, cell.end(), lmbd);
    }

    void all_cell_update() {
        foreach_parallel([](CellPtr& c, size_t i) {
            c->update();
        });
    }
    void update() {
        unsigned int rmsz = remove_index.count();
        auto& rm_arr = remove_index.raw_arr();
        for(unsigned int i=0;i<rmsz;i++){
cell[rm_arr[i]] = cell.back();
cell.pop_back();
        }
        cell.insert(cell.end(), add_cell.raw_arr().begin(), add_cell.raw_arr().begin()+add_cell.count());
        add_cell.clear();
        remove_index.clear();
    }
};

#endif // CELL_H
