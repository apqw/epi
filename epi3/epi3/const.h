#pragma once
#include <vector>
//recursion not allowed
#define OUTPUTDIR "output"
#define THREAD_NUM (8)
template<typename T>
using Arr1D = std::vector<T>;

template<typename T>
using Arr2D = std::vector<std::vector<T>>;

template<typename T>
using Arr3D = std::vector<std::vector<std::vector<T>>>;

static int SYSTEM = 1;
static int WHOLE = 1;
static int BASAL = 0;

static constexpr double EPS = 1e-7;

enum CELL_STATE {
	ALIVE=0,
	DEAD=1,
	DISA=2,
	UNUSED=3, //
	FIX=4,
	BLANK=5,
	DER=6,
	MUSUME=7,
	AIR=8,
	MEMB=9
};
/*
	STATE_NUM:èÛë‘ÇÃéÌóﬁêî
*/
static constexpr int STATE_NUM = 10;
static constexpr int MAX_CELL_CONNECT_NUM = 400;
static constexpr int MAX_CELL_NUM = 30000;
static constexpr double DT_Cell = 0.01;
static constexpr double DT_Ca = 0.02;

/*
LX,LY,LZ:åvéZóÃàÊÇÃÉTÉCÉY
*/
static constexpr double LX = 50.0;
static constexpr double LY = 50.0;
static constexpr double LZ = 100.0;

/*
NX,NY,NZ:åvéZóÃàÊÇÃï™äÑêî
*/
static constexpr int NX = 100;
static constexpr int NY = 100;
static constexpr int NZ = 200;

/*
dx,dy,dz:åvéZóÃàÊÇÃï™äÑïù(L**,N**Ç©ÇÁåvéZÇ≥ÇÍÇÈ)
*/
static constexpr double dx = LX / NX;
static constexpr double dy = LY / NY;
static constexpr double dz = LZ / NZ;

static constexpr double dh = 0.5;

static constexpr double COMPRESS_FACTOR = 4;
static constexpr double P_MEMB = 1 / COMPRESS_FACTOR;
static constexpr double eps_m = 0.01;
static constexpr double DER_DER_CONST = 0.2;
static constexpr double THRESH_SP = 3.0;
static constexpr double K_TOTAL = 3.0;
static constexpr double K_DESMOSOME_RATIO = 0.01;
static constexpr double K_DESMOSOME = K_TOTAL*K_DESMOSOME_RATIO;
static constexpr double LJ_THRESH = 1.2;
static constexpr double Kspring = 25.0;
static constexpr double Kspring_d = 5.0;
static constexpr double R_max = 1.4;
static constexpr double R_der = 1.4;
static constexpr double R_memb = 1.0;
static constexpr double delta_L = 0.01*R_max;
static constexpr double eps_L = 0.14;
static constexpr double delta_R = 0.4*R_der;
static constexpr double para_ljp2 = 0.005;
static constexpr double Kspring_division = 5.0;
static constexpr double agki_max = 6.0;
static constexpr double fac = 1;
static constexpr double agki_max_fix = fac*agki_max;
static constexpr double stoch_div_time_ratio = 0.25;
static constexpr double accel_div = 1.0;
static constexpr double eps_kb = 0.1;
static constexpr double u0 = 0.122;
static constexpr double alpha_b = 5.0;
static constexpr double ADHE_CONST = 31.3;
static constexpr int DISA_conn_num_thresh = 11;
static constexpr double THRESH_DEAD = 22.0;
static constexpr double eps_kk = 0.10*0.5;
static constexpr double eps_ks = 0.10*0.5;
static constexpr double S0 = 0.122*0.2;//u0*0.2;
static constexpr double alpha_k = 2.0;
static constexpr double accel_diff = 1.0;
static constexpr double lipid_rel = 0.05;
static constexpr double ubar = 0.25;
static constexpr double delta_lipid = 0.05;
static constexpr double lipid = 0.032;
static constexpr double unpair_dist_coef = 0.9;
static constexpr double AREA_GRID_ORIGINAL = 2.0;
static constexpr double AREA_GRID = AREA_GRID_ORIGINAL + 1e-7;
static constexpr int ANX = (int)((double)LX / AREA_GRID_ORIGINAL + 0.5);
static constexpr int ANY = (int)((double)LY / AREA_GRID_ORIGINAL + 0.5);
static constexpr int ANZ = (int)((double)LZ / AREA_GRID_ORIGINAL);
static constexpr int N3 = 200; //max grid cell num
static constexpr int N2 = 400; //max conn num
//memb seat size NMX*NMY
static constexpr int NMX = 100;
static constexpr int NMY = 100;
static constexpr int NUM_ITR = (int)1e6;
static constexpr double KBEND = 0.5;

static constexpr int Ca_ITR = 1000;
static constexpr double v0 = 0.97;
static constexpr double p0 = 0;
static constexpr double w0 = 0.99;
static constexpr double a0 = 0;
static constexpr double B0 = 0;

#define KEY_DERMAL_CHANGE 1
#define STOCHASTIC 1

constexpr double constabs(double a);

void static_error_check();

class FObj{};
