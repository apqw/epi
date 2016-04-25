#pragma once
#include <array>
#include <immintrin.h>

template<typename T,unsigned X,unsigned Y,unsigned Z>
using Arr3D = std::array<std::array<std::array<T, Z>, Y>, X>;
template<typename T,unsigned X, unsigned Y>
using Arr2D = std::array<std::array<T, Y>, X>;
template<typename T,unsigned X>
using Arr = std::array<T, X>;

namespace CONST {
	
	constexpr double dt_cell = 0.01;
	constexpr double dt_ca = 0.02;
	
    constexpr int MAX_CELL_NUM = 100000;
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

	//生存時間
	constexpr double THRESH_DEAD = 22.0, DUR_DEAD = 2.0, DUR_ALIVE = 0.5, THRESH_SP = 3.0;

    constexpr double FAC_MAP = 2.0;
    constexpr double R_memb = 1.0;
    constexpr double R_max = 1.4;
    constexpr double R_der=1.4;

	constexpr double KBEND = 0.5;
	constexpr double eps_m =0.01;
	constexpr int N_NEIGHBOR = 4;
	constexpr double COMPRESS_FACTOR= 4;
	constexpr double P_MEMB = 1.0 / COMPRESS_FACTOR;
	constexpr double DER_DER_CONST = 0.2;
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

        constexpr int AIR_STIM=0.1;
	}

}

enum CELL_STATE {
	ALIVE=0,
	DEAD=1,
	DISA=2,
	FIX=4,
	BLANK=5,
	DER=6,
	MUSUME=7,
	AIR=8,
	MEMB=9
};

class ENV {
public:
	int NMEMB;
	int NDER;
	int current_cell_num;
    alignas(32) double Ca2P_value[CONST::NX+1][CONST::NY+1][CONST::NZ+1];
    alignas(32) double B_value[CONST::NX+1][CONST::NY+1][CONST::NZ+1];

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
    alignas(32) int cell_map2[CONST::NZ+1][CONST::NY+1][CONST::NZ+1];
    alignas(32) int cell_radius[CONST::MAX_CELL_NUM];
    alignas(32) double cell_ca_ave[CONST::MAX_CELL_NUM];
    alignas(32) double air_stim[CONST::NX+1][CONST::NY+1][CONST::NZ+1];
    double zzmax;
    alignas(32) int cell_div_pair_index[CONST::MAX_CELL_NUM];
    alignas(32) int origin_stem_cell_index[CONST::MAX_CELL_NUM];//stem cell num?

	alignas(32) int memb_idx_jr[CONST::NMX*CONST::NMY],
		memb_idx_jrr[CONST::NMX*CONST::NMY],
		memb_idx_jl[CONST::NMX*CONST::NMY],
		memb_idx_jll[CONST::NMX*CONST::NMY],
		memb_idx_ju[CONST::NMX*CONST::NMY],
		memb_idx_juu[CONST::NMX*CONST::NMY],
		memb_idx_jb[CONST::NMX*CONST::NMY],
		memb_idx_jbb[CONST::NMX*CONST::NMY];
};
