#pragma once
#include <vector>
#include <array>

#define STOCHASTIC (1)
#define SYSTEM (0)
#define WHOLE (0)
#define BASAL (1)
#define MALIG_NUM (0)
#define OUTPUTDIR "output"
#define CUT (1000)
#define OUT_ENERGY (100)

static const std::string output_dir=OUTPUTDIR;
//#define AGE_DBG
//#define UNPAIR_DBG
template<typename T>
using Arr1D = std::vector<T>;


template<typename T>
using Arr2D = std::vector<Arr1D<T>>;

template<typename T>
using Arr3D = std::vector<Arr2D<T>>;

template<typename T,unsigned Z>
using RawArr1D = std::array<T, Z>;

template<typename T,unsigned Y,unsigned Z>
using RawArr2D = std::array<RawArr1D<T, Z>, Y>;

template<typename T,unsigned X,unsigned Y,unsigned Z>
using RawArr3D = std::array<RawArr2D<T, Y, Z>, X>;


enum CELL_STATE:unsigned int {
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

enum CELL_STATE_MASK:unsigned int {
	ALIVE_M = 1u<<0,
	DEAD_M = 1u << 1,
	DISA_M = 1u << 2,
	UNUSED_M = 1u << 3, //
	FIX_M = 1u << 4,
	BLANK_M = 1u << 5,
	DER_M = 1u << 6,
	MUSUME_M = 1u << 7,
	AIR_M = 1u << 8,
	MEMB_M = 1u << 9
};
/*
inline CELL_STATE_NUM conv_state_num(CELL_STATE state) {
	int targetlevel = 0;
	int index = state;
	while (index >>= 1) ++targetlevel;
	return (CELL_STATE_NUM)targetlevel;
}
*/
namespace cont {
	
	static constexpr int THREAD_NUM = 8;
	static constexpr double ATP_init = 0;
	static constexpr double ext_stim_init = 0;
	static constexpr int STATE_NUM = 9;
	static constexpr int DEFAULT_MAX_CELL_NUM = 30000;

	static constexpr int MAX_CELL_CONNECT_NUM = 400;
	static constexpr int MAX_CELL_NUM = 30000;
	static constexpr double DT_Cell = 0.01;
    static constexpr double DT_Ca = 0.01;//0.02;

	/*
	LX,LY,LZ:ŒvŽZ—Ìˆæ‚ÌƒTƒCƒY
	*/
	static constexpr double LX = 50.0;
	static constexpr double LY = 50.0;
	static constexpr double LZ = 100.0;

	/*
	NX,NY,NZ:ŒvŽZ—Ìˆæ‚Ì•ªŠ„”
	*/
	static constexpr int NX = 100;
	static constexpr int NY = 100;
	static constexpr int NZ = 200;

	/*
	dx,dy,dz:ŒvŽZ—Ìˆæ‚Ì•ªŠ„•(L**,N**‚©‚çŒvŽZ‚³‚ê‚é)
	*/
	static constexpr double dx = LX / NX; static constexpr double inv_dx = NX / LX;
    static constexpr double dy = LY / NY; static constexpr double inv_dy = NY/LY;
    static constexpr double dz = LZ / NZ; static constexpr double inv_dz = NZ/LZ;

    static constexpr double dh = 0.5; //ok
    static constexpr double COMPRESS_FACTOR = 4;//ok
    static constexpr double P_MEMB = 1. / COMPRESS_FACTOR;//ok
    static constexpr double eps_m = 0.01; //ok
    static constexpr double DER_DER_CONST = 0.2; //ok
    static constexpr double THRESH_SP = 3.0; //ok
    static constexpr double K_TOTAL = 3.0;//ok
    static constexpr double K_DESMOSOME_RATIO = 0.01;//ok
    static constexpr double K_DESMOSOME = K_TOTAL*K_DESMOSOME_RATIO;//ok
    static constexpr double LJ_THRESH = 1.2;//ok
    static constexpr double Kspring = 25.0;//ok
    static constexpr double Kspring_d = 5.0;//ok
    static constexpr double R_max = 1.4;//ok
    static constexpr double R_der = 1.4;//ok
    static constexpr double R_memb = 1.0;//ok
    static constexpr double delta_L = 0.01*R_max;//ok

    static constexpr double delta_R = 0.4*R_der;//ok
    static constexpr double para_ljp2 = 0.005;//ok
    static constexpr double Kspring_division = 5.0;//ok
    static constexpr double agki_max = 5.9;//ok
    static constexpr double fac = 1;//ok
    static constexpr double agki_max_fix = fac*agki_max;//ok
    static constexpr double stoch_div_time_ratio = 0.25;//ok
    static constexpr int div_max = 10;//ok
	static constexpr double accel_div = 1.0;
    static constexpr double eps_kb = 0.1;//ok
    static constexpr double u0 = 0.122;//ok
    static constexpr double alpha_b = 5.0;//ok
    static constexpr double ADHE_CONST = 31.3;//ok
    static constexpr int DISA_conn_num_thresh = 11; //Nc ok
    static constexpr double THRESH_DEAD = 22.0;//ok
    static constexpr double eps_kk = 0.10*0.5;//ok
    static constexpr double eps_ks = 0.10*0.5;//ok
    static constexpr double S0 = 0.122*0.2;//u0*0.2; ok
    static constexpr double alpha_k = 2.0;//ok
	static constexpr double accel_diff = 1.0;
    static constexpr double lipid_rel = 0.05;//ok
    static constexpr double ubar = 0.25;//ok
    static constexpr double delta_lipid = 0.05;//ok
    static constexpr double lipid = 0.032;//ok

#ifdef UNPAIR_DBG
	static constexpr double eps_L = 1;
	static constexpr double unpair_dist_coef = 0.9;
#else
    static constexpr double eps_L = 0.14;//ok
	static constexpr double unpair_dist_coef = 0.9;
#endif
    static constexpr double AREA_GRID_ORIGINAL = 2.0;//ok
    static constexpr double AREA_GRID = AREA_GRID_ORIGINAL + 1e-7;//ok
    static constexpr int ANX = (int)((double)LX / AREA_GRID_ORIGINAL + 0.5);//ok
    static constexpr int ANY = (int)((double)LY / AREA_GRID_ORIGINAL + 0.5);//ok
    static constexpr int ANZ = (int)((double)LZ / AREA_GRID_ORIGINAL);//ok
    static constexpr int N3 = 200; //max grid cell num //ok
    static constexpr int N2 = 400; //max conn num //ok
								   //memb seat size NMX*NMY
    static constexpr int NMX = 100;//ok
    static constexpr int NMY = 100;//ok
    static constexpr int NUM_ITR = 2*((int)1e6); //twice
	static constexpr double KBEND = 0.5;

    static constexpr double Ca_avg_time=10.0;
    static int Ca_ITR = (int)(Ca_avg_time/DT_Ca);//Ca_N ok
    static constexpr double v0 = 0.97;//ok
    static constexpr double p0 = 0;//ok
    static constexpr double w0 = 0.99;//ok //0.99->0.1
    static constexpr double a0 = 0;//ok
    static constexpr double B0 = 0;//ok
    static constexpr double ca_0=0.122;//

	static constexpr double Rad = 1.4;
	static constexpr int irx = (int)(Rad*NX / LX);
	static constexpr int iry = (int)(Rad*NY / LY);
	static constexpr int irz = (int)(Rad*NZ / LZ);

    static constexpr double FAC_MAP = 2.0;//ok
    static constexpr double Kpp = 0.3;//ok
    static constexpr double dp = 0.1;//ha? 0.9->0.1 //ok
    static constexpr double kf = 8.1;//ok
    static constexpr double mu0 = 0.567;//ok
    static constexpr double mu1 = 0.1;//ok
    static constexpr double kmu = 0.05;//ok
    static constexpr double para_b = 0.11;//ok
    static constexpr double para_bb = 0.89;//ok
    static constexpr double para_k1 = 0.7;//ok
    static constexpr double gamma = 2.0;//ok
    static constexpr double kg = 0.1;//ok
    static constexpr double beta_zero = 0.02;//ok
    static constexpr double CA_OUT = 1.0;//ok
    static constexpr double beta = CA_OUT*beta_zero;//ok
    static constexpr double kbc = 0.4;//1.0->0.4 //ok
    static constexpr double Hb = 0.01;//ok
    static constexpr double Cout = 1.0;//ha?? 0.4->1.0 //ok
    static constexpr double para_k2 = 0.7;//ok
    static constexpr double thpri = 1.0;//ok
    static constexpr double thgra = 0.2;//ok
    static constexpr double delta_th = 1.0;//ok
    static constexpr double kpa = 4.0;//ok
    static constexpr double Kgra = kpa;//ok
    static constexpr double Kpri = 6.0;//ok
    static constexpr double delta_K = 1.0;//ok
    static constexpr double H0 = 0.5;//ok
    static constexpr double ca2p_du = 0.01; //ha? 0.03->0.01 //ok
    static constexpr double iage_kitei = 0;//ok
    static constexpr double delta_I = 1.5;//ok
    static constexpr double wd = 0.1;//ok
    static constexpr double Da = 1;//ok
    static constexpr double STIM11 = 0.002;//ok
    static constexpr double Kaa = 0.5;//ok
    static constexpr double AIR_STIM = 0.1;//ok
    static constexpr int SW_THRESH = 20;//ok
    static constexpr double DB = 0.0009;//ok
    static constexpr double DUR_ALIVE = 0.5;//ok
    static constexpr double DUR_DEAD = 2.0;//ok
    static constexpr double kb = 0.025;//0.03->0.025
    static constexpr double T_TURNOVER = 6000.0;//
    static constexpr int NUM_SC_INIT = 1;//ok ha/???->19->1
}
