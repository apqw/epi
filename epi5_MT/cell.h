#ifndef CELL_H
#define CELL_H
#include "atomics.h"
#include "define.h"
#include "swapdata.h"
#include <memory>
#include <vector>
#include <unordered_map>

class CellManager;
class Cell;
using CellPtr = std::shared_ptr<Cell>;



class Cell:public std::enable_shared_from_this<Cell>
{

	static SwapData<atomic_double[cont::MAX_CELL_NUM][3]> pos_s;
	static SwapData<double[cont::MAX_CELL_NUM]>ca2p_s;
	static SwapData<double[cont::MAX_CELL_NUM]>IP3_s;
	static SwapData<double[cont::MAX_CELL_NUM]>ex_inert_s;
	static SwapData<double[cont::MAX_CELL_NUM]>agek_s;
	static SwapData<double[cont::MAX_CELL_NUM]>ageb_s;
	static SwapData<double[cont::MAX_CELL_NUM]>ex_fat_s;
	static SwapData<double[cont::MAX_CELL_NUM]>in_fat_s;
	static SwapData<double[cont::MAX_CELL_NUM]>spr_nat_len_s;
	//static SwapData<int[cont::MAX_CELL_NUM]>rest_div_times_s;
	static SwapData<std::unordered_map<Cell*, double>[cont::MAX_CELL_NUM]>gj_s;
	

	size_t index;

	

public:

	static CellManager cells;
	static void pos_swap();
	static void ca2p_swap();
	static void IP3_swap();
	static void ex_inert_swap();
	static void agek_swap();
	static void ageb_swap();
	static void ex_fat_swap();
	static void in_fat_swap();
	static void spr_nat_len_swap();
	static void load_from_file(std::string path);
	static void output(std::string path);

	static unsigned int nmemb;
	static unsigned int nder;
	//static void 
	/*
	 * 位置
	 * 並列で加算など可能
	 * 相互作用の計算で同時に同じ細胞の値に加算する可能性がある
	 */
	SwapArrAccessor2<SwapData<atomic_double[cont::MAX_CELL_NUM][3]>,Cell::pos_s, 0> x;
	SwapArrAccessor2<SwapData<atomic_double[cont::MAX_CELL_NUM][3]>,Cell::pos_s, 1> y;
	SwapArrAccessor2<SwapData<atomic_double[cont::MAX_CELL_NUM][3]>,Cell::pos_s, 2> z;



	/*
	 * 接続された細胞（へのポインタ）のリスト
	 * 追加のみ並列実行できる
	 * タイムステップを越えてポインタを保持しなく、大量に追加することが考えられるので、コストの安い生のポインタを使用
	 */

	Lockfree_push_stack<Cell*, cont::MAX_CONNECT_CELL_NUM> connected_cell;

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
	SwapArrAccessor1<SwapData<double[cont::MAX_CELL_NUM]>, Cell::ca2p_s> ca2p;

	/*
	 * Ca2+濃度(平均)
	 */
	double ca2p_avg;

	/*
	 * IP3
	 */
	SwapArrAccessor1<SwapData<double[cont::MAX_CELL_NUM]>, Cell::IP3_s> IP3;

	/*
	 * 不活性化物質
	 */
	SwapArrAccessor1<SwapData<double[cont::MAX_CELL_NUM]>, Cell::ex_inert_s> ex_inert;

	/*
	 * 分化後の年齢
	 */
	SwapArrAccessor1<SwapData<double[cont::MAX_CELL_NUM]>, Cell::agek_s> agek;

	/*
	 * 分化前の年齢
	 */
	SwapArrAccessor1<SwapData<double[cont::MAX_CELL_NUM]>, Cell::ageb_s> ageb;

	/*
	 * 細胞外脂質
	 */
	SwapArrAccessor1<SwapData<double[cont::MAX_CELL_NUM]>, Cell::ex_fat_s> ex_fat;

	/*
	 * 細胞内脂質
	 */
	SwapArrAccessor1<SwapData<double[cont::MAX_CELL_NUM]>, Cell::in_fat_s> in_fat;

	/*
	 * 細胞分裂中のバネの自然長
	 *
	 * 他の細胞の値も書き換えるが、同時に書き換えないようにできる。
	 */

	SwapArrAccessor1<SwapData<double[cont::MAX_CELL_NUM]>, Cell::spr_nat_len_s> spr_nat_len;

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
	SwapArrAccessor1<SwapData<std::unordered_map<Cell*, double>[cont::MAX_CELL_NUM]>, Cell::gj_s> gj;

	/*
	 * diff_u
	 * Ca2+のみ
	 */
	double diff_u;

	size_t get_index();
	void set_index(size_t i);
	void migrate(size_t dest_idx);

	CELL_STATE state;

	static CellPtr create(CELL_STATE _state, double _x = 0, double _y = 0, double _z = 0,
		double _radius = cont::R_max, double _ca2p = cont::ca2p_init, double _ca2p_avg = cont::ca2p_init,
		double _IP3 = cont::IP3_init, double _ex_inert = cont::ex_inert_init,
		double _agek = 0, double _ageb = 0,
		double _ex_fat = 0, double _in_fat = 0,
		double _spr_nat_len = 0,
		double _div_age_thresh = 0,
		int _rest_div_times = 0,
		bool _is_malignant = false);

	/*
		直接呼ばずにCell::createのほうを使う
	*/
	Cell(CELL_STATE _state,
		double _radius = cont::R_max, double _ca2p_avg = cont::ca2p_init,
		double _div_age_thresh = 0,
		bool _is_malignant = false);

};

#endif // CELL_H
