#ifndef CELL_H
#define CELL_H
#include "atomics.h"
#include "define.h"
#include "swapdata.h"
#include "utils.h"
#include <memory>
#include <vector>
#include <unordered_map>


class Cell:public std::enable_shared_from_this<Cell>
{
	size_t index;
	Cell* _dermis=nullptr;

	struct Memb_data {
		Cell* memb_u; Cell* memb_b; Cell* memb_l; Cell* memb_r;
		Cell* memb_bb; Cell* memb_ll;
		double nv[3]; double ipn;
		double mv[3]; double ipm;
		double dn, dm;
	};

	struct ctor_cookie{};

public:
	friend class CellManager;
	/*
	 * 位置
	 * 並列で加算など可能
	 * 相互作用の計算で同時に同じ細胞の値に加算する可能性がある
	 */
	SwapArrAccessor2<atomic_double[cont::MAX_CELL_NUM][3],0> x;
	SwapArrAccessor2<atomic_double[cont::MAX_CELL_NUM][3],1> y;
	SwapArrAccessor2<atomic_double[cont::MAX_CELL_NUM][3],2> z;



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
	 * 半径
	 * 書き換えないので並列への対応は無し
	 */
	const double radius;

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
	SwapArrAccessor1<double[cont::MAX_CELL_NUM]> ex_inert;

	/*
	 * 分化後の年齢
	 */
	SwapArrAccessor1<double[cont::MAX_CELL_NUM]> agek;

	/*
	 * 分化前の年齢
	 */
	SwapArrAccessor1<double[cont::MAX_CELL_NUM]> ageb;

	/*
	 * 細胞外脂質
	 */
	SwapArrAccessor1<double[cont::MAX_CELL_NUM]> ex_fat;

	/*
	 * 細胞内脂質
	 */
	SwapArrAccessor1<double[cont::MAX_CELL_NUM]> in_fat;

	/*
	 * 細胞分裂中のバネの自然長
	 *
	 * 他の細胞の値も書き換えるが、同時に書き換えないようにできる。
	 */

	SwapArrAccessor1<double[cont::MAX_CELL_NUM]> spr_nat_len;

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

	/*
	 * 格子としての位置
	 */
	int lat[3];

	/*
	 * gj
	 */
	SwapUMapArrAccessor<Cell*,double,cont::MAX_CELL_NUM> gj;

	/*
	 * diff_u
	 * Ca2+のみ
	 */
	double diff_u;

	

	size_t get_index() const;
	void set_index(size_t i);
	
	void set_dermis(Cell* d);
	const Cell* dermis()const;

	void migrate(size_t dest_idx);

	CELL_STATE state;

	double spring_force_to_memb;
	Memb_data md;

	/*
	ctor_cookieがprivateなので自分自身もしくはfriendなclassからのみ生成できる
	*/
	Cell(ctor_cookie,CELL_STATE _state,
		SwapData<atomic_double[cont::MAX_CELL_NUM][3]>& pos_s,
		SwapData<double[cont::MAX_CELL_NUM]>&ca2p_s,
		SwapData<double[cont::MAX_CELL_NUM]>&IP3_s,
		SwapData<double[cont::MAX_CELL_NUM]>&ex_inert_s,
		SwapData<double[cont::MAX_CELL_NUM]>&agek_s,
		SwapData<double[cont::MAX_CELL_NUM]>&ageb_s,
		SwapData<double[cont::MAX_CELL_NUM]>&ex_fat_s,
		SwapData<double[cont::MAX_CELL_NUM]>&in_fat_s,
		SwapData<double[cont::MAX_CELL_NUM]>&spr_nat_len_s,
		SwapData<std::unordered_map<Cell*, double>[cont::MAX_CELL_NUM]>&gj_s,
		double _radius = cont::R_max, double _ca2p_avg = cont::ca2p_init,
		double _div_age_thresh = 0,
		bool _is_malignant = false);

};

#endif // CELL_H
