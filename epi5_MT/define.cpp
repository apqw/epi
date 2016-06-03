#include "define.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <map>
#include <vector>
#include <functional>
#include <algorithm>
#include <cctype>

#define SWITCHER str_switch
#define FORMAT_SWITCHER str_f_switch
#define STORE_SIZE (16)

std::map < std::string, std::function<void(void*)>> SWITCHER;
std::map < std::string, std::string> FORMAT_SWITCHER;
#define pdefd(name) FORMAT_SWITCHER[#name]="%lf";SWITCHER[#name]=[&](void*data){cont::name=*(double*)data;}
#define pdefui(name) FORMAT_SWITCHER[#name]="%u";SWITCHER[#name]=[&](void*data){cont::name=*(unsigned int*)data;}
#define pdefi(name) FORMAT_SWITCHER[#name]="%d";SWITCHER[#name]=[&](void*data){cont::name=*(int*)data;}
namespace cont {
	void set_K_DESMOSOME(double _K_TOTAL, double _K_DESMOSOME_RATIO)
	{
		K_DESMOSOME = _K_TOTAL*_K_DESMOSOME_RATIO;
	}
	void set_delta_L(double _R_max) {
		delta_L = 0.01*_R_max;
	}
	void set_delta_R(double _R_der) {
		delta_R = 0.4*_R_der;
	}
	void set_agki_max_fix(double _fac, double _agki_max) {
		agki_max_fix = _fac*_agki_max;
	}

	void set_Ca_ITR(double _Ca_avg_time, double _DT_Ca) {
		Ca_ITR = (int)(_Ca_avg_time / _DT_Ca);
	}

	void set_beta(double _CA_OUT, double _beta_zero) {
		beta = _CA_OUT*_beta_zero;
	}

	void set_Kgra(double _kpa) {
		Kgra = _kpa;
	}

	void set_ir(double _Rad) {
		irx = (int)(_Rad*NX / LX);
		iry = (int)(_Rad*NY / LY);
		irz = (int)(_Rad*NZ / LZ);
	}

	void set_P_MEMB(double _COMPRESS_FACTOR) {
		P_MEMB = 1. / _COMPRESS_FACTOR;
	}


	void init_param_loader() {
		pdefd(DT_Cell);
		pdefd(DT_Ca);

		pdefd(dh);
		pdefd(COMPRESS_FACTOR);
		//pdefd(P_MEMB);
		pdefd(eps_m);
		pdefd(DER_DER_CONST);
		pdefd(THRESH_SP);
		pdefd(K_TOTAL);
		pdefd(K_DESMOSOME_RATIO);
		
		pdefd(LJ_THRESH);
		pdefd(Kspring);
		pdefd(Kspring_d);
		pdefd(R_max);
		pdefd(R_der);
		pdefd(R_memb);
		pdefd(para_ljp2);
		pdefd(Kspring_division);
		pdefd(agki_max);
		pdefd(fac);
		pdefd(stoch_div_time_ratio);
		pdefui(div_max);
		pdefd(accel_div);
		pdefd(eps_kb);
		pdefd(alpha_b);
		pdefd(ADHE_CONST);
		pdefui(DISA_conn_num_thresh);
		pdefd(eps_kk);
		pdefd(eps_ks);
		pdefd(S0);
		pdefd(alpha_k);
		pdefd(accel_diff);
		pdefd(lipid_rel);
		pdefd(ubar);
		pdefd(delta_lipid);
		pdefd(lipid);
		pdefd(eps_L);
		pdefd(unpair_dist_coef);
		pdefui(N3);
		pdefui(N2);
		pdefui(NMX);
		pdefui(NMY);
		pdefui(NUM_ITR);
		pdefd(KBEND);
		pdefd(Ca_avg_time);
		pdefd(ca2p_init);
		pdefd(u0);
		pdefd(ex_inert_init);
		pdefd(v0);
		pdefd(IP3_init);
		pdefd(v0);
		pdefd(gj_init);
		pdefd(w0);
		pdefd(ATP_init);
		pdefd(a0);
		pdefd(ext_stim_init);
		pdefd(B0);
		pdefd(Rad);
		pdefd(FAC_MAP);
		pdefd(Kpp);
		pdefd(dp);
		pdefd(kf);
		pdefd(mu0);
		pdefd(mu1);
		pdefd(kmu);
		pdefd(para_b);
		pdefd(para_bb);
		pdefd(para_k1);
		pdefd(gamma);
		pdefd(kg);
		pdefd(beta_zero);
		pdefd(CA_OUT);
		pdefd(kbc);
		pdefd(Hb);
		pdefd(Cout);
		pdefd(para_k2);
		pdefd(thpri);
		pdefd(thgra);
		pdefd(delta_th);
		pdefd(kpa);
		pdefd(Kpri);
		pdefd(delta_K);
		pdefd(H0);
		pdefd(ca2p_du);
		pdefd(iage_kitei);
		pdefd(delta_I);
		pdefd(wd);
		pdefd(Da);
		pdefd(STIM11);
		pdefd(Kaa);
		pdefd(AIR_STIM);
		pdefui(AIR_STIM);
		pdefd(DB);
		pdefd(DUR_ALIVE);
		pdefd(DUR_DEAD);
		pdefd(kb);
		pdefd(T_TURNOVER);
		pdefui(NUM_SC_INIT);
		pdefd(stoch_corr_coef);
		pdefi(MALIG_NUM);
	}

	void calc_param()
	{
		set_P_MEMB(COMPRESS_FACTOR);
		set_K_DESMOSOME(K_TOTAL, K_DESMOSOME_RATIO);
		set_delta_L(R_max);
		set_delta_R(R_der);
		set_agki_max_fix(fac, agki_max);
		set_Ca_ITR(Ca_avg_time, DT_Ca);
		set_ir(Rad);
		set_beta(CA_OUT, beta_zero);
		set_Kgra(kpa);
	}

	

	void load_param(std::string path)
	{

		std::ifstream pfstrm;
		pfstrm.open(path, std::ios::in);

		if (!pfstrm) {
			assert(pfstrm);
			throw "load failed";
		}

		std::string line;
		char store[STORE_SIZE]; //sufficiently large size
		while (std::getline(pfstrm, line)) {
			if (line[0] == '#')continue;
			if (line.empty())continue;
			line.erase(std::remove_if(line.begin(), line.end(), std::isspace), line.end());
			std::istringstream stream(line);
			std::string tmp;
			std::vector<std::string> keyval;
			while (std::getline(stream, tmp, '=')) {
				
				keyval.push_back(tmp);
			}
			if (!FORMAT_SWITCHER.count(keyval[0])) {
				printf("unrecognized parameter:%s\n", keyval[0].c_str());
				continue;
			}
			for (int i = 0; i < STORE_SIZE; i++) {
				store[i] = 0;
			}
			if (sscanf(keyval[1].c_str(),FORMAT_SWITCHER[keyval[0]].c_str(), store) != 1) {
				printf("parsing value of %s failed. value=%s\n", keyval[0].c_str(),keyval[1].c_str());
				continue;
			}
			SWITCHER[keyval[0]](store);
		}
		calc_param();
		
		printf("param initialized %lf %u\n",cont::KBEND,NMX);
	}
}