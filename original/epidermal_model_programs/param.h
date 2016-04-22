#define K_TOTAL 3.0
#define K_DESMOSOME_RATIO 0.01

#define NAME_ENERGY_OUTPUT "energy.dat"
#define SW_THRESH 20
#define INIT_SW (SW_THRESH-1) // SW_THRESH
#define NUM_ITR 1e6
#define Nc 11
#define CUT 1000
#define OUT_ENERGY 100
#define NX 100 // 100 // grid for calculation of heat equations
#define NY 100
#define NZ 200
#define DBG 0
#define DBGP 1000
#define DBG_CD 0
#define DBG_LJ 0

#define div_max 10 

#define N_NEIGHBOR 4 
#define NN_NEIGHBOR 4
#define OUTPUTDIR "output"


//#define OUTPUTDIR "/media/disk/kobayashi/basal/output"

#define NN 30000  // max #cells //done
#define N2 400 // used in lj[1..NN][1..N2], wj　//done
#define N3 200 // area

#define LENGTH_X 50.0
#define LENGTH_Y 50.0
#define LENGTH_Z 100.0

#define eps_m 0.01

#define DT 0.01 // for cell dynamics
#define Dtca 0.02 // for calcium dynamics

#define para_b 0.11 //do
#define para_bb 0.89
#define beta_zero 0.02
#define gamma 2.
#define para_k1 0.7 //do
#define para_k2 0.7 //do
#define para_kg 0.1 //do
#define para_kf 8.1 //done
#define para_kp 0.2
#define para_kmu 0.05 //do
#define mu0 0.567 //do
#define mu1 0.1 //do
#define CA_OUT 1.0 //do
#define beta (CA_OUT*beta_zero) //do

#define AMAX 0.002//1.
#define AMIN 0.
#define UMAX 0.5 
#define UMIN 0.1 
#define VMAX 1.
#define VMIN 0.
#define PMAX 0.5
#define PMIN 0.

#define Rad 1.4     /* 0.02 */
#define a0 0.
#define para_delta 0.1
#define lez2 (2. * Rad)
#define p0 0.
#define u0 0.122
#define S0 0.122 * 0.2
#define v0 0.97
#define w0 0.99//0.9999
#define pc 0.

//#define LJ_THRESH 1.2 // threshold for LJ connection
extern double LJ_THRESH;
extern int SYSTEM;
#define LJ_EPSILON 0.01

#define FAC_MAP 2.0

#define para_kau 0.02
#define para_I0 1.
#define ubar 0.25
#define para_kaa 0.5//0.2
#define para_Kgra para_kpa //done
#define para_Kpri 6. //para_kpa * 10. //done
#define para_H0 0.5 //done
#define para_kpp 0.3 //done

#define eps_w 0.001

#define delta_K 1.0//0.1  //0.1 //done
#define delta_I 1.5//0.15 //0.01 //done

#define para_iage_kitei 0. //done

#define T_TURNOVER 0.0
#define KEY_FORCED_SC 0
#define para_alpha_b 5.0
#define para_alpha_k 2.0 

#define PERIODIC 1

#define R_memb 1.0
#define R_max 1.4 
#define R_der 1.4 

#define LX LENGTH_X //done
#define LY LENGTH_Y //done
#define LZ LENGTH_Z //done

#define AREA_GRID_ORIGINAL 2.0 // grid size for area, default:2.0
#define AREA_GRID (AREA_GRID_ORIGINAL+1e-7) // grid size for area, default:2.0
#define ANX (int)((double)LX / AREA_GRID_ORIGINAL+0.5)
#define ANY (int)((double)LY / AREA_GRID_ORIGINAL+0.5)
#define ANZ (int)((double)LZ / AREA_GRID_ORIGINAL)

#define para_agki_max 6.0 // 3.0
#define fac 1//1
#define para_agki_max_fix (fac*para_agki_max) // 3*para_agki_max 

#define stoch_div_time_ratio 0.25 /* ratio of stochastic div. time average to deterministic div. time */

#define ubar_re 0.40

#define para_lipid 0.032 //0.0039
#define para_lipid_rel 0.05 //0.15

#define THRESH_DEAD 22.0 //done
#define DUR_DEAD 2.0	//done
#define DUR_ALIVE 0.5	//done
#define SIGMA_DEAD 0.05
#define SIGMA_ALIVE 0.01

#define STIM11 0.002 //done
	
#define THRESH_SP 3.0 //7.0 // 9.0 // age of stratum spinosum ( a.k.a. T_pri) //done

#define para_k 1. 
#define para_kx 1.
#define para_ky 1.
#define para_kz 1.0  

#define para_mu 1. 
#define para_er 0.05

#define para_rc R
#define para_rc2 R2
#define para_test 10.


/*saibougai*/
/* state list */
#define ALIVE 0 /* 生存 */
#define DEAD  1 /* 死亡(角層) */
#define DISA  2 /* 消滅 */
#define FIX   4 /* 基底 */
#define BLANK 5 /* 空白 */ 
#define DER 6   /* 真皮 */
#define MUSUME 7/* 娘幹細胞 */
#define AIR   8 /* "air cells */
#define MEMB 9  /*basal membrane */

extern int NDER, NMEMB, NTRI;
extern int KEY_DERMAL_CHANGE;

extern double dx, dy, dz;
#define dh 0.5

extern int ca_N ;// 1000;
extern double para_thgra ;// 0.2; 
extern double para_thpri ;// 1.0;
extern double delta_th ;// 1.0;
extern double FMAX;//2.0;
extern double FMIN;//0.0;
extern double lez1 ;// (2. * Rad);
extern double lez ;// 5 * Rad;
extern double lbx ;// LX * 0.5;
extern double p0c ;// 0.;
//extern double ue ;// 1.5;//1.5;

extern double Da ;// 1.; //done
extern double du ;// 0.03; //done
extern double dp ;// 0.9; //done
extern double para_kpa ;// 4.;
extern double para_wd ;// 0.1;//0.2;
extern double delta_cos ;// 0.1;

extern char str_para[256] ;// {};

extern double para_dpri ;// 0.;

//extern double eps_kb ;// 0.10 ;
extern double eps_kb ;// 0.10 ;
extern double eps_ks ;// 0.10 * 0.5;      //細胞の状態に応じて状態変数の変動率を変える？
extern double eps_kk ;// 0.10 * 0.5;

extern double B0 ;// 0.0;

extern double DB ;// 0.0009;    //Bの計算で用いる変化率
extern double delta_IR ;// 3.;
extern double para_kb ;// 0.03; //done

extern double para_kbc ;// 1.0;
extern double para_Hb ;// 0.01;

extern double para_Cout ;// 0.4;

/*-------------------*/
extern double para_Kfat ;// 0.05;
extern double calory ;// 2.0;

extern double myu ;// 0.05;

extern double GF;//0.003;

extern int irx, iry, irz;  

/*******/

/*角層が剥がれる状態変数値*/
#define ADHE_CONST 31.3

extern double debug_double;
extern int debug_int;

extern double AIR_STIM;

#define eps_L 0.14
#define delta_L (0.01*R_max)

#define para_ljp2 0.005
#define delta_R (0.4*R_der)

extern double accel_div, accel_diff;
extern int MALIGNANT;

extern double DER_DER_CONST;
extern double DD_FAC;
extern double DD_LJ_FAC;
extern double PRESS_CONST;
extern double KBEND;
extern double distmax, distmin;

#define COMPRESS_FACTOR 4
#define P_MEMB (1.0/COMPRESS_FACTOR) // compression ratio for the membrane

#define NMX 100
#define NMY 100
// NMX * NMY must agree with NMEMB
#define NUM_SC_INIT 1

#define WHOLE 0
#define BASAL 1
