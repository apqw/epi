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


class Cell:public std::enable_shared_from_this<Cell>
{
	static constexpr double agki_max		= 6.0;
	static constexpr double fac				= 1;
	static constexpr double agki_max_fix	= fac*agki_max;
	static constexpr int	div_max			= 20;


	static double get_div_age_thresh(CELL_STATE state);
	static int correct_div_times(CELL_STATE state, int given_times);

	size_t index;
	const Cell* _dermis=nullptr;

	struct Memb_data {
		Cell* memb_u; Cell* memb_b; Cell* memb_l; Cell* memb_r;
		Cell* memb_bb; Cell* memb_ll;
		double nv[3]; double ipn;
		double mv[3]; double ipm;
		double dn, dm;
	};

	struct ctor_cookie{};

    bool removed=false;
public:
	friend class CellManager;

        CELL_STATE state;

        const int fix_origin;
	/*
	 * 位置
	 * 並列で加算など可能
	 * 相互作用の計算で同時に同じ細胞の値に加算する可能性がある
	 */
	dual_double x;
	dual_double y;
	dual_double z;



	/*
	 * 接続された細胞（へのポインタ）のリスト
	 * 追加のみ並列実行できる
	 * タイムステップを越えてポインタを保持しなく、大量に追加することが考えられるので、コストの安い生のポインタを使用
	 */

	Lockfree_push_stack<Cell*, cont::MAX_CONNECT_CELL_NUM> connected_cell;

	/*
	 * 分裂中のペアへのポインタ
	 */
	CellPtr pair=nullptr;



	/*
	 * Ca2+濃度(平均化されていない)
	 * この値はCa2+の計算部分以外では使われないはず
	 */
	SwapArrAccessor1<double[cont::MAX_CELL_NUM]> ca2p;

	/*
	 * Ca2+濃度(平均)
	 */
	double ca2p_avg;

	/*
	 * IP3
	 */
	SwapArrAccessor1<double[cont::MAX_CELL_NUM]> IP3;

	/*
	 * 不活性化物質
	 */
	double ex_inert;

	/*
	 * 分化後の年齢
	 */
	double agek;

	/*
	 * 分化前の年齢
	 */
	double ageb;

	/*
	 * 細胞外脂質
	 */
	double ex_fat;

	/*
	 * 細胞内脂質
	 */
	double in_fat;

	/*
	 * 細胞分裂中のバネの自然長
	 *
	 * 他の細胞の値も書き換えるが、同時に書き換えないようにできる。
	 */

	double spr_nat_len;

    /*
     * 半径
     * 書き換える
     */
    double radius;

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
	int rest_div_times;

	/*
	 * 悪性かどうか
	 */
	bool is_malignant;

    bool is_touch;

	/*
	 * 格子としての位置
	 */
	int lat[3];

	/*
	 * gj
	 */
    std::unordered_map<const Cell*,gj_value> gj;

	/*
	 * diff_u
	 * Ca2+のみ
	 */
	double diff_u;




	inline size_t get_index() const
	{
		return index;
	}
	void set_index(size_t i);
	
	void set_dermis(const Cell* d);
	const Cell* dermis()const;

	void migrate(size_t dest_idx);



	double spring_force_to_memb;
	Memb_data md;

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

inline double p_cell_dist_sq(const Cell*const RESTRICT c1, const Cell*const RESTRICT c2)
{
	return p_dist_sq(c1->x(), c1->y(), c1->z(), c2->x(), c2->y(), c2->z());
}

inline bool no_double_count(const Cell*const c1,const Cell*const c2) {
	return c1->get_index() > c2->get_index();
}

#endif // CELL_H
