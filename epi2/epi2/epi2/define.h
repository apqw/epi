#pragma once
#include <array>
template<typename T,unsigned X,unsigned Y,unsigned Z>
using Arr3D = std::array<std::array<std::array<T, Z>, Y>, X>;
template<typename T,unsigned X, unsigned Y>
using Arr2D = std::array<std::array<T, Y>, X>;
template<typename T,unsigned X>
using Arr = std::array<T, X>;

namespace CONST {
	
	constexpr double dt_cell = 0.01;
	constexpr double dt_ca = 0.02;
	
	constexpr int MAX_CELL_NUM = 3000;
	constexpr int MAX_REACT_CELL_NUM = 400;

	//計算領域のサイズ
	constexpr double Lx = 50.0, Ly = 50.0, Lz = 100.0;

	//計算領域の分割数
	constexpr int NX = 100, NY = 100, NZ = 200;

	//グリッドサイズ
	constexpr double dx = Lx / NX, dy = Ly / NY, dz = Lz / NZ;

	//生存時間
	constexpr double THRESH_DEAD = 22.0, DUR_DEAD = 2.0, DUR_ALIVE = 0.5, THRESH_SP = 3.0;

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
	Arr3D<double, CONST::NX, CONST::NY, CONST::NZ> Ca2P_value;
	Arr3D<double, CONST::NX, CONST::NY, CONST::NZ> B_value;

	Arr<double, CONST::MAX_CELL_NUM> cell_x;
	Arr<double, CONST::MAX_CELL_NUM> cell_y;
	Arr<double, CONST::MAX_CELL_NUM> cell_z;
	Arr<int, CONST::MAX_CELL_NUM> cell_x_index;
	Arr<int, CONST::MAX_CELL_NUM> cell_y_index;
	Arr<int, CONST::MAX_CELL_NUM> cell_z_index;
	Arr<double,CONST::MAX_CELL_NUM> cell_connected_num;
	Arr2D<double, CONST::MAX_CELL_NUM, CONST::MAX_REACT_CELL_NUM> cell_connected_index; //max==cell_connected_num
	Arr<double, CONST::MAX_CELL_NUM> cell_P;
	Arr<double, CONST::MAX_CELL_NUM> cell_c;
	Arr<double, CONST::MAX_CELL_NUM> cell_h;
	Arr<double, CONST::MAX_CELL_NUM> cell_ageb;
	Arr<double, CONST::MAX_CELL_NUM> cell_agek;
	Arr<CELL_STATE, CONST::MAX_CELL_NUM> cell_state;
	Arr2D<double,CONST::MAX_CELL_NUM, CONST::MAX_REACT_CELL_NUM> cell_w;

};
