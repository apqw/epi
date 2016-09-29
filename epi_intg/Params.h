/*
 * Params.h
 *
 *  Created on: 2016/09/20
 *      Author: yasu7890v
 */
#pragma once
#ifndef PARAMS_H_
#define PARAMS_H_
#include "define.h"
#include "utils.h"
#include <vector>
#include <type_traits>
#include <typeinfo>
#include <initializer_list>
#include <map>
#include <stdexcept>
#include <string>
#include <array>
static const constexpr std::array<const char*, 2> str_sign = { "\"","\'" };

template<typename T>
T str_cast(const std::string& str) {
    return T();
}

template<>
unsigned int str_cast<unsigned int>(const std::string& str);

template<>
int str_cast<int>(const std::string& str);

template<>
double str_cast<double>(const std::string& str);

template<>
float str_cast<float>(const std::string& str);

template<>
bool str_cast<bool>(const std::string& str);

template<>
std::string str_cast<std::string>(const std::string& str);

template<typename VarTy>
VarTy try_get_with_cast(const std::string& item, const char* first = "NONE") {
    try {
        return str_cast<VarTy>(item);
    }
    catch (std::invalid_argument& ie) {
        throw std::invalid_argument("Failed to read the value of "_s + first + "(='" + item + "') with the following exception:\n" + typeid(ie).name() + ": " + ie.what() + "\n");
    }
    catch (std::out_of_range& oe) {
        throw std::out_of_range("Failed to read the value of "_s + first + "(='" + item + "') with the following exception:\n" + typeid(oe).name() + ": " + oe.what() + "\n");
    }
    catch (std::exception& e) {
        throw std::logic_error("Failed to read the value of "_s + first + "(='" + item + "') with the following exception:\n" + typeid(e).name() + ": " + e.what() + "\n");
    }
}

template<typename VarTy, typename K, typename V, typename K2>
VarTy force_vget(const std::map<K, V>& mp, const std::vector<K2>& kv) {
    for (auto& key : kv) {
        auto it = mp.find(key);
        if (it != mp.end()) {
            return try_get_with_cast<VarTy>(it->second, key);
        }
    }
    throw std::logic_error(std::string(kv[0]) + " must be specified.");
}

struct ParamInfo {
private:
    enum TypeEnum {
        UINT, INT, REAL,BOOL,STRING
    };


    void* __param;
    TypeEnum pTy;
    template<typename Ty, typename Nm>
    ParamInfo(Ty& _p, std::initializer_list<Nm>_n, typename std::enable_if<std::is_same<Ty, unsigned int>::value>::type* = nullptr) :__param(&_p), name(_n), pTy(TypeEnum::UINT) {
    }

    template<typename Ty, typename Nm>
    ParamInfo(Ty& _p, std::initializer_list<Nm>_n, typename std::enable_if<std::is_same<Ty, int>::value>::type* = nullptr) : __param(&_p), name(_n), pTy(TypeEnum::INT) {
    }
    template<typename Ty, typename Nm>
    ParamInfo(Ty& _p, std::initializer_list<Nm>_n, typename std::enable_if<std::is_same<Ty, real>::value>::type* = nullptr) : __param(&_p), name(_n), pTy(TypeEnum::REAL) {
    }
    template<typename Ty, typename Nm>
    ParamInfo(Ty& _p, std::initializer_list<Nm>_n, typename std::enable_if<std::is_same<Ty, bool>::value>::type* = nullptr) : __param(&_p), name(_n), pTy(TypeEnum::BOOL) {
    }

    template<typename Ty, typename Nm>
    ParamInfo(Ty& _p, std::initializer_list<Nm>_n, typename std::enable_if<std::is_same<Ty, std::string>::value>::type* = nullptr) : __param(&_p), name(_n), pTy(TypeEnum::STRING) {
    }

public:
    std::vector<const char*> name;

    
    
    template<typename Ty, typename... Nm,class = typename std::enable_if<sizeof...(Nm)>=1>::type>
    ParamInfo(Ty& _p, Nm... _n) :ParamInfo(_p, { _n... }) {}

    template<typename K, typename V>
    void set_from(const std::map<K, V>& mp) {
        bool chk = true;
        switch (pTy) {
        case UINT:

            *reinterpret_cast<unsigned int*>(__param) = force_vget<unsigned int>(mp, name);
            break;
        case INT:
            *reinterpret_cast<int*>(__param) = force_vget<int>(mp, name);
            break;
        case REAL:
            *reinterpret_cast<real*>(__param) = force_vget<real>(mp, name);
            break;
        case BOOL:
            *reinterpret_cast<bool*>(__param) = force_vget<bool>(mp, name);
            break;
        case STRING:
            *reinterpret_cast<std::string*>(__param) = force_vget<std::string>(mp, name);
            break;
        default:
            break;
        }
    }

    std::string get_as_string()const {
        bool chk = true;
        switch (pTy) {
        case UINT:
            return std::to_string(*reinterpret_cast<unsigned int*>(__param));
            break;
        case INT:
            return std::to_string(*reinterpret_cast<int*>(__param));
            break;
        case REAL:
            return std::to_string(*reinterpret_cast<real*>(__param));
            break;
        case BOOL:
            return *reinterpret_cast<bool*>(__param) ? "1" : "0";
            break;
        case STRING:
            return std::string(str_sign[0])+*reinterpret_cast<std::string*>(__param)+str_sign[0];
            break;
        default:
            return "!ERROR";
            break;
        }
    }

};

class Params {
private:
    std::vector<ParamInfo> piset;

    void calculate_param();
public:
	//no const for params plz
	//use const for this class

	//general,misc,uncategorized
	 unsigned int MEMB_NUM_X;
	 unsigned int MEMB_NUM_Y;

	 unsigned int MAX_CONNECT_CELL_NUM; //use for warning
	// unsigned int MAX_CELL_NUM;
	 real DT_Cell;
	 real DT_Ca;
	 unsigned int MALIG_NUM;
	 real LX,LY,LZ;
	 unsigned int NX,NY,NZ;
	 real R_max,R_der,R_memb;
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
	 unsigned int ANX,ANY,ANZ;
	 unsigned int N2,N3;//for warning

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

     unsigned int SYSTEM;
     unsigned int CUT;
     bool NEW_BEND_POT;

     std::string outputdir;
     
	 Params();
     //Params(const char* paramfile);
     Params(const std::string& paramfile);
     void generate_paramfile(const char* out_path)const;
	 void load(const std::string& paramfile);
};


#endif /* PARAMS_H_ */
