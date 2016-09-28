#ifndef DEFINE_H
#define DEFINE_H
//#define NDEBUG
#include <string>
#include <memory>
#define cdefd static constexpr double
#define cdefi static constexpr int
#define cdefui static constexpr unsigned int
#define cdefs static constexpr size_t

#define defd static double
#define rdefd static double&
#define defi static int
#define rdefi static int&
#define defui static unsigned int
#define rdefui static unsigned int&
#define defs static size_t
#define rdefs static size_t&

#define cmax(a,b) ((a)>(b)?(a):(b))
#ifdef _WIN32
#define RESTRICT __restrict
#else
#define RESTRICT __restrict__
#endif

class Cell;
using CellPtr = Cell*;//std::shared_ptr<Cell>;

using cmask_ty = uint_fast8_t;
using CELL_STATE_internal = uint_fast8_t;

/*
	DO NOT CHANGE 
*/
enum CELL_STATE:CELL_STATE_internal {
    ALIVE = 0,
    DEAD = 1,
    DISA = 2,
    UNUSED = 3, //
    FIX = 4,
    BLANK = 5,
    DER = 6,
    MUSUME = 7,
    AIR = 8,
    MEMB = 9
};

enum CELL_STATE_MASK :uint_fast16_t {
	ALIVE_M = 1u<<ALIVE,
	DEAD_M = 1u << DEAD,
	DISA_M = 1u<<DISA,
	UNUSED_M = 1u<<UNUSED, //
	FIX_M = 1u<<FIX,
	BLANK_M = 1u<<BLANK,
	DER_M = 1u<<DER,
	MUSUME_M = 1u<<MUSUME,
	AIR_M = 1u<<AIR,
	MEMB_M = 1u<<MEMB
};

enum BoundaryType {
	DIRICHLET, NEUMANN, PERIODIC
};

#define TRI_MEMB 1
//#define DIAG_BEND 1
//no diag
namespace cont {
	
    //NMX
    static constexpr unsigned int NMX = 100;//ok
    //NMY
    static constexpr unsigned int NMY = 100;//ok

    //細胞の状態の数
	cdefs STATE_NUM = 10;

    //細胞の最大接続数
	cdefs MAX_CONNECT_CELL_NUM = 400;

    //細胞の最大数
    cdefs MAX_CELL_NUM = NMX*NMY+50000*10;

    //細胞の時間ステップ
	cdefd DT_Cell = 0.01;

    //カルシウムの時間ステップ
	cdefd DT_Ca = 0.01;//0.02;

    //これ以下のインデックスの幹細胞が悪性扱い
	cdefi MALIG_NUM = 0;

    /*
        LX,LY,LZ:計算領域のサイズ
    */
    cdefd LX = 50.0;
    cdefd LY = 50.0;
	cdefd LZ = 100.0;

	/*
	NX,NY,NZ:計算領域の分割数
	*/
    cdefs NX = 100;
    cdefs NY = 100;
	cdefs NZ = 200;

	/*
	dx,dy,dz:計算領域の分割幅(L**,N**から計算される)
	*/
	cdefd dx = LX / NX; cdefd inv_dx = NX / LX;
	cdefd dy = LY / NY; cdefd inv_dy = NY / LY;
	cdefd dz = LZ / NZ; cdefd inv_dz = NZ / LZ;

    //細胞の(最大)半径
    cdefd R_max = 1.4;//ok

    //derの半径
    cdefd R_der = 1.4;//ok

    //基底膜の細胞の半径
    cdefd R_memb = 1.0;//ok

    //基底膜の細胞の間隔の圧縮比
    cdefd COMPRESS_FACTOR = 4;//ok
   // cdefui MEMB_ADHE_RANGE = (unsigned int)cmax((int)((R_max/R_der)*COMPRESS_FACTOR/2)-1,0);//test

    //幹細胞と基底膜を接着する際の接着範囲(0で最も近いもののみ)
    cdefui MEMB_ADHE_RANGE=0;

    //THRESH_SP
    cdefd THRESH_SP = 3.0; //ok 3.0->1.0
	
    //LJ_THRESH
	cdefd LJ_THRESH = 1.2;//ok

    //基底膜の細胞間の接続数(現在は上下左右の4つ)
    /*
    cdefui MEMB_ADJ_CONN_NUM=
        #ifdef DIAG_BEND
            8
        #else
            4
        #endif
            ;
    */
    cdefui MEMB_ADJ_CONN_NUM=
        #ifdef TRI_MEMB
            6
        #else
            4
        #endif
            ;

    //THRESH_DEAD
	cdefd THRESH_DEAD = 22.0;//ok

    //全体のループ回数
    cdefui NUM_ITR = 4 * ((int)1e6);

    //Ca2+の計算(平均化)に使う時間
	cdefd Ca_avg_time = 10.0;
	
    //Ca2+濃度の初期値
	cdefd ca2p_init = 0.122;

    //不活性化物質濃度の初期値
	cdefd ex_inert_init = 0.97;

    //IP3濃度の初期値
	cdefd IP3_init = 0;

    //GJ発現率の初期値
	cdefd gj_init = 0.99;

    //ATPの初期値
	cdefd ATP_init = 0;

    //外刺激物質濃度の初期値
	cdefd ext_stim_init = 0;

	
    //FAC_MAP
	cdefd FAC_MAP = 2.0;//ok
	
    //SW_THRESH
	cdefi SW_THRESH = 20;//ok

    //T_TURNOVER
    cdefd T_TURNOVER = 6000.0;

    //NUM_SC_INIT
    cdefi NUM_SC_INIT = 1;

	cdefd stoch_corr_coef = 1;
    
}
struct cell_stat{
    unsigned long long k_cornified_timestep;
    unsigned long long k_disap_timestep;
    unsigned long long k_aging_start_timestep;
};

//SYSTEM
cdefui SYSTEM = 0;

//WHOLE
cdefui WHOLE = 0;

//BASAL
cdefui BASAL = 1;

//CUT
cdefui CUT=1000;


static constexpr auto OUTPUTDIR="output";

static constexpr auto last_data_uvp_name="last_data_uvp";
static constexpr auto last_data_cell_name="last_data_cell";
static constexpr auto last_data_w_name="last_data_w_alt";
static constexpr auto last_data_a_name="last_data_a_alt";
static constexpr auto last_data_B_name="last_data_B_alt";
static constexpr auto stat_data_name="cell_stat";

//確率的に分裂するか
static constexpr bool STOCHASTIC = true;

//新しい曲げ弾性エネルギーを使うかどうか
#define NEW_BEND_POT 1
//static constexpr bool FORCE_CORNIF = true;


static constexpr double __sqrt3 =1.73205080757;

static constexpr double y_comp_ratio=__sqrt3*0.5;

static constexpr size_t NMY_tri = (((size_t)(cont::NMY/y_comp_ratio))/(size_t)2)*2;

template<typename T, size_t X, size_t Y, size_t Z>
class Field;
template<typename T>
using FArr3D = Field<T, cont::NX + 1, cont::NY + 1, cont::NZ + 1>;

template<typename T>
using RawArr3D = T[cont::NX + 1][cont::NY + 1][cont::NZ + 1];
#endif // DEFINE_H
