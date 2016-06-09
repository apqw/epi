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

#ifdef _WIN32
#define RESTRICT __restrict
#else
#define RESTRICT __restrict__
#endif

class Cell;
using CellPtr = Cell*;//std::shared_ptr<Cell>;



/*
	DO NOT CHANGE 
*/
enum CELL_STATE:uint_fast8_t {
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

namespace cont {
	
	cdefs STATE_NUM = 10;
	cdefs MAX_CONNECT_CELL_NUM = 400;
	cdefs MAX_CELL_NUM = 30000;
	cdefd DT_Cell = 0.01;
	cdefd DT_Ca = 0.01;//0.02;
	cdefi MALIG_NUM = 0;
										 /*
										 LX,LY,LZ:ŒvŽZ—Ìˆæ‚ÌƒTƒCƒY
										 */
	cdefd LX = 50.0;
	cdefd LY = 50.0;
	cdefd LZ = 100.0;

	/*
	NX,NY,NZ:ŒvŽZ—Ìˆæ‚Ì•ªŠ„”
	*/
	cdefs NX = 100;
	cdefs NY = 100;
	cdefs NZ = 200;

	/*
	dx,dy,dz:ŒvŽZ—Ìˆæ‚Ì•ªŠ„•(L**,N**‚©‚çŒvŽZ‚³‚ê‚é)
	*/
	cdefd dx = LX / NX; cdefd inv_dx = NX / LX;
	cdefd dy = LY / NY; cdefd inv_dy = NY / LY;
	cdefd dz = LZ / NZ; cdefd inv_dz = NZ / LZ;

	cdefd dh = 0.5; //ok
	cdefd COMPRESS_FACTOR = 4;//ok
	//void set_P_MEMB(double _COMPRESS_FACTOR);
	cdefd P_MEMB = 1. / COMPRESS_FACTOR;//ok
	cdefd eps_m = 0.01; //ok
	cdefd DER_DER_CONST = 0.2; //ok
	cdefd THRESH_SP = 3.0; //ok 3.0->1.0
	cdefd K_TOTAL = 3.0;//ok
	cdefd K_DESMOSOME_RATIO = 0.01;//ok
	cdefd K_DESMOSOME = K_TOTAL*K_DESMOSOME_RATIO;//ok
	//void set_K_DESMOSOME(double _K_TOTAL,double _K_DESMOSOME_RATIO);
	cdefd LJ_THRESH = 1.2;//ok
	cdefd Kspring = 25.0;//ok
	cdefd Kspring_d = 5.0;//ok
	cdefd R_max = 1.4;//ok
	cdefd R_der = 1.4;//ok
	cdefd R_memb = 1.0;//ok
	cdefd delta_L = 0.01*R_max;//ok
	//void set_delta_L(double _R_max);

	cdefd delta_R= 0.4*R_der;//ok
	//void set_delta_R(double _R_der);
	cdefd para_ljp2 = 0.005;//ok
	cdefd Kspring_division = 5.0;//ok
	cdefd agki_max = 6.0;//ok
	cdefd fac = 1;//ok
	cdefd agki_max_fix = fac*agki_max;//ok
	//void set_agki_max_fix(double _fac, double _agki_max);
	cdefd stoch_div_time_ratio = 0.25;//ok
	cdefui div_max = 15;//ok
	cdefd accel_div = 1.0;
	cdefd eps_kb = 0.1;//ok

	
	cdefd alpha_b = 5.0;//ok
	cdefd ADHE_CONST = 31.3;//ok
	cdefui DISA_conn_num_thresh = 11; //Nc ok
	//rdefui Nc = DISA_conn_num_thresh;
	cdefd THRESH_DEAD = 22.0;//ok
	cdefd eps_kk = 0.10*0.5;//ok
	cdefd eps_ks = 0.10*0.5;//ok
	cdefd S0 = 0.122*0.2;//u0*0.2; ok
	cdefd alpha_k = 2.0;//ok
	cdefd accel_diff = 1.0;
	cdefd lipid_rel = 0.05*2.0;//ok
	cdefd ubar = 0.25;//ok
	cdefd delta_lipid = 0.05;//ok
	cdefd lipid = 0.032*2.0;//ok

#ifdef UNPAIR_DBG
	cdefd eps_L = 1;
	cdefd unpair_dist_coef = 0.9;
#else
	cdefd eps_L = 0.14;//ok
	cdefd unpair_dist_coef = 0.9;
#endif
	cdefd AREA_GRID_ORIGINAL = 2.0;//ok
	cdefd AREA_GRID = AREA_GRID_ORIGINAL + 1e-7;//ok
	cdefs ANX = (size_t)((double)LX / AREA_GRID_ORIGINAL + 0.5);//ok
	cdefs ANY = (size_t)((double)LY / AREA_GRID_ORIGINAL + 0.5);//ok
	cdefs ANZ = (size_t)((double)LZ / AREA_GRID_ORIGINAL);//ok
	cdefui N3 = 200; //max grid cell num //ok
	cdefui N2 = 400; //max conn num //ok
								   //memb seat size NMX*NMY
	cdefui NMX = 100;//ok
	cdefui NMY = 100;//ok
	cdefui NUM_ITR = 2 * ((int)1e6); //twice
	cdefd KBEND = 0.5;

	cdefd Ca_avg_time = 10.0;
	cdefi Ca_ITR = (int)(Ca_avg_time / DT_Ca);//Ca_N ok
	//void set_Ca_ITR(double _Ca_avg_time, double _DT_Ca);

	cdefd ca2p_init = 0.122;
	//rdefd u0 = ca2p_init;//ok

	cdefd ex_inert_init = 0.97;
	//rdefd v0 = ex_inert_init;//ok

	cdefd IP3_init = 0;
	//rdefd p0 = IP3_init;//ok

	cdefd gj_init = 0.99;
	//rdefd w0 = gj_init;//ok //0.99->0.1

	cdefd ATP_init = 0;
	//rdefd a0 = ATP_init;//ok

	cdefd ext_stim_init = 0;
	//rdefd B0 = ext_stim_init;//ok

	cdefd Rad = 1.4;
	//void set_ir(double _Rad);
	
	cdefi irx = (int)(Rad*NX / LX);
	cdefi iry= (int)(Rad*NY / LY);
	cdefi irz = (int)(Rad*NZ / LZ);

	cdefd FAC_MAP = 2.0;//ok
	cdefd Kpp = 0.3;//ok
	cdefd dp = 0.1;//ha? 0.9->0.1 //ok
	cdefd kf = 8.1;//ok
	cdefd mu0 = 0.567;//ok
	cdefd mu1 = 0.1;//ok
	cdefd kmu = 0.05;//ok
	cdefd para_b = 0.11;//ok
	cdefd para_bb = 0.89;//ok
	cdefd para_k1 = 0.7;//ok
	cdefd gamma = 2.0;//ok
	cdefd kg = 0.1;//ok
	cdefd beta_zero = 0.02;//ok
	cdefd CA_OUT = 1.0;//ok
    cdefd beta = CA_OUT*beta_zero;//ok
	//void set_beta(double _CA_OUT, double _beta_zero);
	cdefd kbc = 0.4*1.2;//1.0->0.4 //ok //0.6
	cdefd Hb = 0.01;//ok
	cdefd Cout = 1.0;//ha?? 0.4->1.0 //ok
	cdefd para_k2 = 0.7;//ok
	cdefd thpri = 1.0;//ok
	cdefd thgra = 0.2;//ok
	cdefd delta_th = 1.0;//ok
	cdefd kpa = 4.0;//ok
	cdefd Kgra = kpa;//ok
	//void set_Kgra(double _kpa);
	cdefd Kpri = 6.0;//ok
	cdefd delta_K = 1.0;//ok
	cdefd H0 = 0.5;//ok
	cdefd ca2p_du = 0.01; //ha? 0.03->0.01 //ok
	cdefd iage_kitei = 0;//ok
	cdefd delta_I = 1.5;//ok
	cdefd wd = 0.1;//ok
	cdefd Da = 1;//ok
	cdefd STIM11 = 0.002;//ok
	cdefd Kaa = 0.5;//ok
	cdefd AIR_STIM = 0.1;//ok
	cdefi SW_THRESH = 20;//ok
	cdefd DB = 0.0009;//ok
	cdefd DUR_ALIVE = 0.5;//ok
	cdefd DUR_DEAD = 2.0;//ok
	cdefd kb = 0.025;//0.03->0.025
	cdefd T_TURNOVER = 6000.0;//
	cdefi NUM_SC_INIT = 1;//ok ha/???->19->1
	cdefd stoch_corr_coef = 1;
}

cdefui SYSTEM = 0;
cdefui WHOLE = 0;
cdefui BASAL = 1;
static constexpr bool STOCHASTIC = true;

static constexpr bool FORCE_CORNIF = true;

template<typename T, size_t X, size_t Y, size_t Z>
class Field;
template<typename T>
using FArr3D = Field<T, cont::NX + 1, cont::NY + 1, cont::NZ + 1>;

template<typename T>
using RawArr3D = T[cont::NX + 1][cont::NY + 1][cont::NZ + 1];
#endif // DEFINE_H
