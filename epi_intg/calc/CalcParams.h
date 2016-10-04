#pragma once
#include "../Params.h"

class CalcParams:public Params {
private:
   // std::vector<ParamInfo> piset;

    void calculate_param();
public:
    //no const for params plz
    //use const for this class

    //general,misc,uncategorized
    typedef Params Base;
    unsigned int MEMB_NUM_X;
    unsigned int MEMB_NUM_Y;

    unsigned int MAX_CONNECT_CELL_NUM; //use for warning
                                       // unsigned int MAX_CELL_NUM;
    real DT_Cell;
    real DT_Ca;
    unsigned int MALIG_NUM;
    real LX, LY, LZ;
    unsigned int NX, NY, NZ;
    real R_max, R_der, R_memb;
    real COMPRESS_FACTOR;
    unsigned int MEMB_ADHE_RANGE;

    real THRESH_SP;
    real LJ_THRESH;
    unsigned int MEMB_ADJ_CONN_NUM;
    real THRESH_DEAD;
    unsigned int NUM_ITR;

    real Ca_avg_time;
    real ca2p_init;
    real ex_inert_init;
    real IP3_init;
    real gj_init;
    real ATP_init;
    real ext_stim_init;

    real FAC_MAP;
    unsigned int SW_THRESH;
    real T_TURNOVER;
    unsigned int NUM_SC_INIT;

    bool USE_TRI_MEMB;

    //Cell
    real agki_max;
    real fac;
    real agki_max_fix;//fac*agki_max
    unsigned int div_max;

    //Ca2+
    unsigned int Ca_ITR;
    real Kpp;
    real dp;

    real kbc;
    real Hb;
    real Cout;

    real kg;
    real gamma;

    real mu0;
    real mu1;
    real kmu;

    real para_b;
    real para_bb;
    real para_k1;

    real k_flux;
    real beta_zero;
    real CA_OUT;
    real leak_from_storage;

    real thgra;
    real delta_th;
    real thpri;

    real kpa;
    real Kgra;
    real delta_K;
    real Kpri;

    real delta_I;
    real iage_kitei;

    real para_k2;

    real H0;

    real wd;

    real ca2p_du;

    real STIM11;
    real Kaa;

    real Da;
    real AIR_STIM;

    //connection
    real AREA_GRID_ORIGINAL;
    real AREA_GRID;
    unsigned int ANX, ANY, ANZ;
    unsigned int N2, N3;//for warning

                        //interaction
    real eps_m;

    real P_MEMB;
    real DER_DER_CONST;
    real K_TOTAL;
    real K_DESMOSOME_RATIO;
    real K_DESMOSOME;
    real Kspring;
    real Kspring_d;

    real KBEND;

    real delta_R; real delta_R_coef;

    real para_ljp2;
    real Kspring_division;

    //state renewal
    real stoch_div_time_ratio;
    real S2;

    real accel_div;
    real eps_kb;

    real alpha_b;

    real S0;

    real eps_ks;
    real accel_diff;

    real alpha_k;

    real eps_kk;
    real S1;

    real ubar;
    real delta_lipid;
    real delta_sig_r1;
    real lipid_rel;
    real lipid;

    real delta_L; real delta_L_coef;

    real ADHE_CONST;
    real DISA_conn_num_thresh;

    real eps_L;
    real unpair_dist_coef;

    //ext stim
    real kb;
    real DUR_ALIVE;
    real DUR_DEAD;
    real DB;

    real y_tri_comp_ratio;
    unsigned int SYSTEM;
    unsigned int CUT;
    bool NEW_BEND_POT;

    std::string outputdir;
    
    CalcParams();
    //Params(const char* paramfile);
    CalcParams(const std::string& paramfile);
    //void generate_paramfile(const char* out_path)const;
    void init();

    void load(const std::string& paramfile);
    
};
