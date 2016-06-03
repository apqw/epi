#ifndef DEFINE_H
#define DEFINE_H
#include <string>

#define cdefd static constexpr double
#define cdefi static constexpr int
#define cdefui static constexpr unsigned int
#define cdefs static constexpr size_t

#define defd static double
#define rdefd static double&
#define defi static int
#define defui static unsigned int
#define defs static size_t

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

enum BoundaryType {
	DIRICHLET, NEUMANN, PERIODIC
};


/*
namespace cont{
cdefs MAX_CONNECTED_CELL_NUM=400;
defs MAX_CELL_NUM=30000;
cdefd R_max = 1.4;//ok
cdefd R_der = 1.4;//ok
cdefd R_memb = 1.0;//ok
cdefd ca2p_init=0.122;
cdefd IP3_init=0;
cdefd ex_inert_init=0.97;

cdefd agki_max = 6.0;//ok
cdefd fac = 1;//ok
cdefd agki_max_fix = fac*agki_max;//ok

cdefui div_max = 10;
cdefui MALIG_NUM = 0;

cdefui NX = 100;
cdefui NY = 100;
cdefui NZ = 200;
}
*/
namespace cont {
	void load_param(std::string path);
	void init_param_loader();
	void calc_param();
	

	cdefs MAX_CONNECT_CELL_NUM = 400;
	cdefs MAX_CELL_NUM = 30000;
	defd DT_Cell = 0.01;
	defd DT_Ca = 0.01;//0.02;
	defi MALIG_NUM = 0;
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

	defd dh = 0.5; //ok
	defd COMPRESS_FACTOR = 4;//ok
	void set_P_MEMB(double _COMPRESS_FACTOR);
	defd P_MEMB;// = 1. / COMPRESS_FACTOR;//ok
	defd eps_m = 0.01; //ok
	defd DER_DER_CONST = 0.2; //ok
	defd THRESH_SP = 3.0; //ok 3.0->1.0
	defd K_TOTAL = 3.0;//ok
	defd K_DESMOSOME_RATIO = 0.01;//ok
	defd K_DESMOSOME;// = K_TOTAL*K_DESMOSOME_RATIO;//ok
	void set_K_DESMOSOME(double _K_TOTAL,double _K_DESMOSOME_RATIO);
	defd LJ_THRESH = 1.2;//ok
	defd Kspring = 25.0;//ok
	defd Kspring_d = 5.0;//ok
	defd R_max = 1.4;//ok
	defd R_der = 1.4;//ok
	defd R_memb = 1.0;//ok
	defd delta_L;// = 0.01*R_max;//ok
	void set_delta_L(double _R_max);

	defd delta_R;// = 0.4*R_der;//ok
	void set_delta_R(double _R_der);
	defd para_ljp2 = 0.005;//ok
	defd Kspring_division = 5.0;//ok
	defd agki_max = 6.0;//ok
	defd fac = 1;//ok
	defd agki_max_fix;// = fac*agki_max;//ok
	void set_agki_max_fix(double _fac, double _agki_max);
	defd stoch_div_time_ratio = 0.25;//ok
	defui div_max = 15;//ok
	defd accel_div = 1.0;
	defd eps_kb = 0.1;//ok

	
	defd alpha_b = 5.0;//ok
	defd ADHE_CONST = 31.3;//ok
	defui DISA_conn_num_thresh = 11; //Nc ok
	defd THRESH_DEAD = 22.0;//ok
	defd eps_kk = 0.10*0.5;//ok
	defd eps_ks = 0.10*0.5;//ok
	defd S0 = 0.122*0.2;//u0*0.2; ok
	defd alpha_k = 2.0;//ok
	defd accel_diff = 1.0;
	defd lipid_rel = 0.05*2.0;//ok
	defd ubar = 0.25;//ok
	defd delta_lipid = 0.05;//ok
	defd lipid = 0.032*2.0;//ok

#ifdef UNPAIR_DBG
	defd eps_L = 1;
	defd unpair_dist_coef = 0.9;
#else
	defd eps_L = 0.14;//ok
	defd unpair_dist_coef = 0.9;
#endif
	cdefd AREA_GRID_ORIGINAL = 2.0;//ok
	cdefd AREA_GRID = AREA_GRID_ORIGINAL + 1e-7;//ok
	cdefs ANX = (size_t)((double)LX / AREA_GRID_ORIGINAL + 0.5);//ok
	cdefs ANY = (size_t)((double)LY / AREA_GRID_ORIGINAL + 0.5);//ok
	cdefs ANZ = (size_t)((double)LZ / AREA_GRID_ORIGINAL);//ok
	defui N3 = 200; //max grid cell num //ok
	defui N2 = 400; //max conn num //ok
								   //memb seat size NMX*NMY
	defui NMX = 100;//ok
	defui NMY = 100;//ok
	defui NUM_ITR = 2 * ((int)1e6); //twice
	defd KBEND = 0.5;

	defd Ca_avg_time = 10.0;
	defi Ca_ITR;// = (int)(Ca_avg_time / DT_Ca);//Ca_N ok
	void set_Ca_ITR(double _Ca_avg_time, double _DT_Ca);

	defd ca2p_init = 0.122;
	rdefd u0 = ca2p_init;//ok

	defd ex_inert_init = 0.97;
	rdefd v0 = ex_inert_init;//ok

	defd IP3_init = 0;
	rdefd p0 = IP3_init;//ok

	defd gj_init = 0.99;
	rdefd w0 = gj_init;//ok //0.99->0.1

	defd ATP_init = 0;
	rdefd a0 = ATP_init;//ok

	defd ext_stim_init = 0;
	rdefd B0 = ext_stim_init;//ok

	defd Rad = 1.4;
	void set_ir(double _Rad);
	
	defi irx;// = (int)(Rad*NX / LX);
	defi iry;// = (int)(Rad*NY / LY);
	defi irz;// = (int)(Rad*NZ / LZ);

	defd FAC_MAP = 2.0;//ok
	defd Kpp = 0.3;//ok
	defd dp = 0.1;//ha? 0.9->0.1 //ok
	defd kf = 8.1;//ok
	defd mu0 = 0.567;//ok
	defd mu1 = 0.1;//ok
	defd kmu = 0.05;//ok
	defd para_b = 0.11;//ok
	defd para_bb = 0.89;//ok
	defd para_k1 = 0.7;//ok
	defd gamma = 2.0;//ok
	defd kg = 0.1;//ok
	defd beta_zero = 0.02;//ok
	defd CA_OUT = 1.0;//ok
	defd beta;// = CA_OUT*beta_zero;//ok
	void set_beta(double _CA_OUT, double _beta_zero);
	defd kbc = 0.4*1.2;//1.0->0.4 //ok //0.6
	defd Hb = 0.01;//ok
	defd Cout = 1.0;//ha?? 0.4->1.0 //ok
	defd para_k2 = 0.7;//ok
	defd thpri = 1.0;//ok
	defd thgra = 0.2;//ok
	defd delta_th = 1.0;//ok
	defd kpa = 4.0;//ok
	defd Kgra;// = kpa;//ok
	void set_Kgra(double _kpa);
	defd Kpri = 6.0;//ok
	defd delta_K = 1.0;//ok
	defd H0 = 0.5;//ok
	defd ca2p_du = 0.01; //ha? 0.03->0.01 //ok
	defd iage_kitei = 0;//ok
	defd delta_I = 1.5;//ok
	defd wd = 0.1;//ok
	defd Da = 1;//ok
	defd STIM11 = 0.002;//ok
	defd Kaa = 0.5;//ok
	defd AIR_STIM = 0.1;//ok
	defi SW_THRESH = 20;//ok
	defd DB = 0.0009;//ok
	defd DUR_ALIVE = 0.5;//ok
	defd DUR_DEAD = 2.0;//ok
	defd kb = 0.025;//0.03->0.025
	defd T_TURNOVER = 6000.0;//
	defi NUM_SC_INIT = 1;//ok ha/???->19->1
	defd stoch_corr_coef = 1;
}

cdefui SYSTEM = 0;
cdefui WHOLE = 0;
cdefui BASAL = 1;

#endif // DEFINE_H
