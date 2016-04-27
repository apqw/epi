#pragma once
#define MEXP 607
#include <array>
#include <immintrin.h>
#define isNOTMINUS(val) ((((~(val))&0x80000000)>>31) & 1)
#define isNOTZERO(n) ((((n) | ((~(n)) + 1)) >> 31) & 1)
#define WHOLE (0)
#define BASAL (1)
template<typename T,unsigned X,unsigned Y,unsigned Z>
using Arr3D = std::array<std::array<std::array<T, Z>, Y>, X>;
template<typename T,unsigned X, unsigned Y>
using Arr2D = std::array<std::array<T, Y>, X>;
template<typename T,unsigned X>
using Arr = std::array<T, X>;
static int c_SYSTEM = 0;
enum CELL_STATE {
	CELL_UNDEF = 0,
	ALIVE = (1u << 0),
	DEAD = (1u << 1),
	DISA = (1u << 2),
	FIX = (1u << 3),
	BLANK = (1u << 4),
	DER = (1u << 5),
	MUSUME = (1u << 6),
	AIR = (1u << 7),
	MEMB = (1u << 8)
};

namespace CONST {
	constexpr int NUM_ITR = (int)(1E+6);
	constexpr unsigned int ALIVE_FIX_MUSUME = ALIVE | FIX | MUSUME;
	constexpr double dt_cell = 0.01;
	constexpr double dt_ca = 0.02;
	
    constexpr int MAX_CELL_NUM = 30000;
	constexpr int MAX_REACT_CELL_NUM = 400;

	//計算領域のサイズ
	constexpr double Lx = 50.0, Ly = 50.0, Lz = 100.0;

	//計算領域の分割数
	constexpr int NX = 100, NY = 100, NZ = 200;
	constexpr int NMX = 100, NMY = 100;

	

	//グリッドサイズ

	constexpr double dx = Lx / NX, dy = Ly / NY, dz = Lz / NZ;
    constexpr double dxSq=dx*dx,dySq=dy*dy,dzSq=dz*dz;
    constexpr double inv_dxSq = 1/dxSq,inv_dySq=1/dySq,inv_dzSq=1/dzSq;
    constexpr double inv_dx = 1/dx,inv_dy=1/dy,inv_dz=1/dz;

	constexpr double AREA_GRID_ORIGINAL = 2.0;
	constexpr double AREA_GRID = AREA_GRID_ORIGINAL + 1e-7;//?
	constexpr double AREA_GRID_inv = 1 / AREA_GRID;
	constexpr int ANX = (int)((double)Lx / AREA_GRID_ORIGINAL + 0.5);
	constexpr int ANY = (int)((double)Ly / AREA_GRID_ORIGINAL + 0.5);
	constexpr int ANZ = (int)((double)Lx / AREA_GRID_ORIGINAL);
	constexpr int cell_on_lat_max = 200;

	constexpr int N2 = 400;

	//生存時間
	constexpr double THRESH_DEAD = 22.0, DUR_DEAD = 2.0, DUR_ALIVE = 0.5, THRESH_SP = 3.0;

    constexpr double FAC_MAP = 2.0;
    constexpr double R_memb = 1.0;
    constexpr double R_max = 1.4;
    constexpr double R_der=1.4;
	constexpr double delta_R = 0.4*R_der;

	constexpr double KBEND = 0.5;
	constexpr double eps_m =0.01;
	constexpr int N_NEIGHBOR = 4;
	constexpr double COMPRESS_FACTOR= 4;
	constexpr double P_MEMB = 1.0 / COMPRESS_FACTOR;
	constexpr double DER_DER_CONST = 0.2;
	constexpr double ljp2 = 0.005;
	constexpr double K_TOTAL = 3.0;
	constexpr double K_DESMOSOME_RATIO = 0.01;
	constexpr double K_DESMOSOME = K_TOTAL*K_DESMOSOME_RATIO;
	constexpr double LJ_THRESH = 1.2;
	constexpr double Kspring = 25.0;
	constexpr double Kspring_d = 5.0;
	constexpr double Kspring_division = 5;
	constexpr int KEY_DERMAL_CHANGE = 1;
	constexpr double agki_max = 6.0; constexpr double agki_max_inv = 1 / agki_max;
	constexpr double fac = 1.0;
	constexpr double agki_max_fix = fac*agki_max; constexpr double agki_max_fix_inv = 1 / agki_max_fix;
	constexpr double stoch_div_time_ratio = 0.25; constexpr double stoch_div_time_ratio_inv = 1/stoch_div_time_ratio;

	constexpr double accel_div = 1.0;
	constexpr double eps_kb = 0.10;
	constexpr double u0 = 0.122;
	constexpr int MALIGNANT_NUM = 0;
	constexpr int STOCHASTIC = 1;

	constexpr double delta_L = 0.01*R_max;
	constexpr int div_max = 10;
	constexpr double alpha_b = 5.0;
	constexpr double ADHE_CONST = 31.3;
	constexpr int Nc = 11;
	constexpr double lipid_rel = 0.05;
	constexpr double ubar = 0.25;
	constexpr double delta_lipid = 0.05; constexpr double delta_lipid_inv = 1 / delta_lipid;
	constexpr double eps_L = 0.14;
	constexpr double div_fin_coef = 0.9;
	namespace Ca2P {
		//拡散係数
		constexpr double dA = 1.0, dP = 0.1, dc = 0.01, dB = 0.0009;

		constexpr double Kaa = 0.5, Kpp = 0.3, Kbb = 0.025;

		constexpr double Kac = 0.002, Kf = 8.1, Kmu = 0.05, K1 = 0.7, Kg = 0.1, Kbc = 0.4;

		constexpr double mu0 = 0.567, mu1 = 0.1, alpha0 = 0.11;

		constexpr double sub1_alpha0 = 1.0 - alpha0;
		constexpr double gamma = 2.0;
		constexpr double beta_zero = 0.02;
		constexpr double CA_OUT = 1.0;
		constexpr double beta = beta_zero*CA_OUT;

		constexpr double Hb = 0.01;
		constexpr double H0 = 0.5;
		constexpr double Cout = 1.0;

		constexpr double K2 = 0.7; constexpr double K2Sq = K2*K2;
		constexpr double wd0 = 0.1;
		constexpr double eps_w0 = 0.1;

		constexpr double Sk1 = THRESH_DEAD - DUR_ALIVE, Sk2 = THRESH_DEAD + DUR_DEAD, Ss = THRESH_SP;
		constexpr double tau_g = 0.2, tau_s = 1.0;
		constexpr double delta_tau = 1.0, delta_I = 1.5, delta_k = 1.0;
		constexpr double kg = 4.0, ks = 6.0;
		
		constexpr double _r = 1.0;

		constexpr double iage_kitei = 0.0;
		
		//smooth step
		constexpr int ca_N = 1000;

        constexpr double AIR_STIM=0.1;
	}

}



class ENV {
public:
	int NMEMB;
	int NDER;
	int current_cell_num;
	int pair_count;
	int born_num;
	int vacant_num;
	int next_pair_index;
	int disap_num; //dead
	int sw;//???
	alignas(32) int vanished_cell_index[CONST::MAX_CELL_NUM]; //indexed by vacant_num (1 start)
	alignas(32) double tmp_diffu_map[CONST::NX + 1][CONST::NY + 1][CONST::NZ + 1];
    alignas(32) double Ca2P_value[CONST::NX+1][CONST::NY+1][CONST::NZ+1];
    alignas(32) double B_value[CONST::NX+1][CONST::NY+1][CONST::NZ+1];
	alignas(32) int cell_div_times[CONST::MAX_CELL_NUM];
	alignas(32) double cell_fat[CONST::MAX_CELL_NUM];
	alignas(32) double cell_Vc[CONST::MAX_CELL_NUM];
	alignas(32) bool cell_generated[CONST::MAX_CELL_NUM];

    alignas(32) double cell_x[CONST::MAX_CELL_NUM];
    alignas(32) double cell_y[CONST::MAX_CELL_NUM];
    alignas(32) double cell_z[CONST::MAX_CELL_NUM];
    alignas(32) int cell_x_index[CONST::MAX_CELL_NUM];
    alignas(32) int cell_y_index[CONST::MAX_CELL_NUM];
    alignas(32) int cell_z_index[CONST::MAX_CELL_NUM];
    alignas(32) int cell_connected_num[CONST::MAX_CELL_NUM];
    alignas(32) int cell_connected_index[CONST::MAX_CELL_NUM][CONST::MAX_REACT_CELL_NUM]; //max==cell_connected_num
    alignas(32) double cell_P[CONST::MAX_CELL_NUM];
    alignas(32) double cell_c[CONST::MAX_CELL_NUM];
    alignas(32) double cell_h[CONST::MAX_CELL_NUM];
    alignas(32) double cell_ageb[CONST::MAX_CELL_NUM];
    alignas(32) double cell_agek[CONST::MAX_CELL_NUM];
    alignas(32) CELL_STATE cell_state[CONST::MAX_CELL_NUM];
    alignas(32) double cell_w[CONST::MAX_CELL_NUM][CONST::MAX_REACT_CELL_NUM];
    alignas(32) double cell_diff_c[CONST::MAX_CELL_NUM];
    alignas(32) int cell_map[CONST::NX+1][CONST::NY+1][CONST::NZ+1];
    alignas(32) double cell_map2[CONST::NZ+1][CONST::NY+1][CONST::NZ+1];
    alignas(32) int cell_radius[CONST::MAX_CELL_NUM];
    alignas(32) double cell_ca_ave[CONST::MAX_CELL_NUM];
    alignas(32) double air_stim[CONST::NX+1][CONST::NY+1][CONST::NZ+1];
    double zzmax;
    alignas(32) int cell_div_pair_index[CONST::MAX_CELL_NUM];
    alignas(32) int origin_stem_cell_index[CONST::MAX_CELL_NUM];//stem cell num?
	alignas(32) int cell_pair_lesser[CONST::MAX_CELL_NUM];
	alignas(32) int cell_pair_larger[CONST::MAX_CELL_NUM];
	alignas(32) double cell_spring_len[CONST::MAX_CELL_NUM];
	alignas(32) int memb_idx_jr[CONST::NMX*CONST::NMY],
		memb_idx_jrr[CONST::NMX*CONST::NMY],
		memb_idx_jl[CONST::NMX*CONST::NMY],
		memb_idx_jll[CONST::NMX*CONST::NMY],
		memb_idx_ju[CONST::NMX*CONST::NMY],
		memb_idx_juu[CONST::NMX*CONST::NMY],
		memb_idx_jb[CONST::NMX*CONST::NMY],
		memb_idx_jbb[CONST::NMX*CONST::NMY];
};
