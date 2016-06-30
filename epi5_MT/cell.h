#ifndef CELL_H
#define CELL_H
#include "atomics.h"
#include "define.h"
#include "swapdata.h"
#include "utils.h"
#include <memory>
#include <vector>
#include <unordered_map>
#include "cell_conn_value.h"
#include "DualValue.h"


/**
 *  細胞を記述するクラス
 *
 *  細胞のステートごとに派生クラスを作成しようかと考えたが、
 *  ステートの変更が面倒になるのでしていない。
 */
class Cell:public std::enable_shared_from_this<Cell>
{
    /** MUSUME用分裂開始年齢のしきい値*/
	static constexpr double agki_max		= 6.0;

    /** FIX用分裂開始年齢の倍率 */
	static constexpr double fac				= 1;

    /** FIX用分裂開始年齢のしきい値 */
	static constexpr double agki_max_fix	= fac*agki_max;

    /** 最大分裂回数 */
    static constexpr int	div_max			= 15;


	static double get_div_age_thresh(CELL_STATE state);
	static int correct_div_times(CELL_STATE state, int given_times);

	size_t index;
	const Cell* _dermis=nullptr;

    /** 基底膜としてのデータ */
	struct Memb_data {
        /** 奥のmemb */
		Cell* memb_u;
        /** 手前のmemb */
        Cell* memb_b;
        /** 左のmemb */
        Cell* memb_l;
        /** 右のmemb */
        Cell* memb_r;
        /** 2つ奥のmemb */
		Cell* memb_bb; 
        /** 2つ左のmemb */
        Cell* memb_ll;

		double nv[3]; double ipn;
		double mv[3]; double ipm;
		double dn, dm;

        double nv_a[3]; double ipn_a;
        double mv_a[3]; double ipm_a;
        double dn_a, dm_a;
        Cell** adj_memb[4]={&memb_u,&memb_l,&memb_b,&memb_r};

	};


    template<unsigned init_group,unsigned as,unsigned rest_depth,class...Args,typename =typename std::enable_if<sizeof...(Args)==0 && rest_depth >= 1>::type>
    inline bool _memb_exist_sec(const Cell* const target,Args...)const{
        if(target==this)return true;
        //if(current_depth>=max_depth)return false;
        if((*(md.adj_memb[as]))->_memb_exist_sec<init_group,as,rest_depth-1>(target)){
            return true;
        }
        if(init_group==as){
            constexpr int another=(as+1)%4;
            if((*(md.adj_memb[another]))->_memb_exist_sec<init_group,another,rest_depth-1>(target)){
                return true;
            }
        }
        return false;
    }

    template<unsigned init_group,unsigned as,unsigned rest_depth,class...Args,typename =typename std::enable_if<sizeof...(Args)==0 && rest_depth <= 0,void >::type>
    inline bool _memb_exist_sec(const Cell* const target,Args*...)const{
        return target==this;
    }

    /** cookie struct.
     コンストラクタの呼び出し元を制限するのに使う */
	struct ctor_cookie{};
    void migrate(size_t dest_idx);
    void swap_index(Cell* dest);
public:
	friend class CellManager;

    /** 細胞のステート */
    CELL_STATE state;

    /** 分裂元の幹細胞のインデックス */
    const int fix_origin;
	/**
	 * x座標。
	 * 相互作用の計算で同時に同じ細胞の値を操作する可能性がある。
	 */
	dual_double x;

    /**
     * y座標。
     * 並列で加算など可能。
     * 相互作用の計算で同時に同じ細胞の値を操作する可能性がある。
     */
	dual_double y;

    /**
     * z座標。
     * 並列で加算など可能。
     * 相互作用の計算で同時に同じ細胞の値を操作する可能性がある。
     */
	dual_double z;


	/**
	 * 接続された細胞（へのポインタ）のリスト。
	 * 追加のみ並列実行できる。
	 */

	Lockfree_push_stack<Cell*, cont::MAX_CONNECT_CELL_NUM> connected_cell;

	/**
	 * 分裂中のペアへのポインタ。
	 */
	Cell* pair=nullptr;


	/**
	 * Ca2+濃度(平均化されていない)。
	 * この値はCa2+の計算部分以外では使われないはず。
	 */
	SwapArrAccessor1<double[cont::MAX_CELL_NUM]> ca2p;

	/**
	 * Ca2+濃度(平均)
	 */
	double ca2p_avg;

	/**
	 * IP3濃度
	 */
	SwapArrAccessor1<double[cont::MAX_CELL_NUM]> IP3;

	/**
	 * 不活性化物質
	 */
	double ex_inert;

	/**
	 * 分化年齢
	 */
	double agek;

	/**
	 * 分裂年齢
	 */
	double ageb;

	/**
	 * 細胞外脂質
	 */
	double ex_fat;

	/**
	 * 細胞内脂質
	 */
	double in_fat;

	/**
	 * 細胞分裂中のバネの自然長
	 *
	 * 他の細胞の値も書き換えるが、同時に書き換えないようにできる。
	 */

	double spr_nat_len;

    /**
     * 半径
     */
    const double radius;

	/**
	 * 分裂開始年齢のしきい値
	 *
	 * この値に確率に関する係数がかかって分裂開始年齢になる
	 */
	const double div_age_thresh;

	/**
	 * 残りの分裂回数
	 */
	int rest_div_times;

	/**
	 * 悪性かどうか
	 */
	bool is_malignant;

    /** 局在化可能細胞かどうか (データ出力時に初めて設定される)*/
    bool is_touch;

	/**
	 * 格子としての位置
	 */
	int lat[3];

	/**
	 * gj
	 */
    std::unordered_map<const Cell*,gj_value> gj;

	/**
	 * diff_u
     *
	 * Ca2+処理内でのみ使用される
	 */
	double diff_u;

    cell_stat cst;
    /** インデックスを取得 */
	inline size_t get_index() const
	{
		return index;
	}
	void set_index(size_t i);
	
	void set_dermis(const Cell* d);
	const Cell* dermis()const;

	



	double spring_force_to_memb;
	Memb_data md;







    template<unsigned max_depth,typename...Args,typename=typename std::enable_if<sizeof...(Args)==0&&max_depth>=1>::type>
    bool memb_exist(const Cell* const target,Args...)const{
        if(target==this)return true;

#define MFIND(n) if(_memb_exist_sec<n,n,max_depth-1>(target))return true
        MFIND(0);MFIND(1);MFIND(2);MFIND(3);
#undef MFIND
        return false;
    }

    template<unsigned max_depth,typename...Args,typename=typename std::enable_if<sizeof...(Args)==0&&max_depth<=0>::type>
    bool memb_exist(const Cell* const target,Args*...)const{
        return target==this;
    }
	/*
	  ctor_cookieがprivateなので自分自身もしくはfriendなclassからのみ生成できる
	*/
    Cell(ctor_cookie,CELL_STATE _state,int fix_origin,
		SwapData<double[cont::MAX_CELL_NUM]>&ca2p_s,
		SwapData<double[cont::MAX_CELL_NUM]>&IP3_s,
		double _ex_inert=0,
		double _agek = 0,double _ageb = 0,double _ex_fat=0,double _in_fat=0,double _spr_nat_len=0,
		double _x=0, double _y = 0, double _z = 0,int _rest_div_time=0,
		double _radius = cont::R_max, double _ca2p_avg = cont::ca2p_init,
		bool _is_malignant = false);

    void* operator new(size_t s);
    void operator delete(void* p);

};

/** 周期境界条件を考慮した細胞間距離 */
inline double p_cell_dist_sq(const Cell*const RESTRICT c1, const Cell*const RESTRICT c2)
{
	return p_dist_sq(c1->x(), c1->y(), c1->z(), c2->x(), c2->y(), c2->z());
}

/** ダブルカウントの回避条件 */
inline bool no_double_count(const Cell*const c1,const Cell*const c2) {
	return c1->get_index() > c2->get_index();
}

#endif // CELL_H
