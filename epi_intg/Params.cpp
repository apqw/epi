/*
 * Params.cpp
 *
 *  Created on: 2016/09/20
 *      Author: yasu7890v
 */

#include "Params.h"
#include "parser.h"
#include <iostream>

#include <fstream>


#include <cassert>
#include <regex>


template<>
unsigned int str_cast<unsigned int>(const std::string& str) {
    return std::stoul(str);
}

template<>
int str_cast<int>(const std::string& str) {
    return std::stoi(str);
}

template<>
double str_cast<double>(const std::string& str) {
    return std::stod(str);
}

template<>
float str_cast<float>(const std::string& str) {
    return std::stof(str);
}

template<>
bool str_cast<bool>(const std::string& str) {
    return std::stoi(str) == 1;
}

template<>
std::string str_cast<std::string>(const std::string& str) {

    for (const char* s : str_sign) {
        std::regex re("^"_s + s + "(.*)" + s + "$");
        std::smatch match;
        if (std::regex_match(str, match, re)) {
            return match[1];
        }
    }
    throw std::logic_error("The following string data:" + str + " is ill-formatted.");
}

#define gpa(name,...) ParamInfo(name,#name,__VA_ARGS__)
#define gp1(name) ParamInfo(name,#name)

Params::Params():piset{
    gpa(MEMB_NUM_X,"NMX"),
    gpa(MEMB_NUM_Y,"NMY"),
    gp1(MAX_CONNECT_CELL_NUM),
    gp1(DT_Cell),
    gp1(DT_Ca),
    gpa(MALIG_NUM,"MALIGNANT"),
    gp1(LX),gp1(LY),gp1(LZ),
    gp1(NX),gp1(NY),gp1(NZ),
    gp1(R_max),gp1(R_der),gp1(R_memb),
    gp1(COMPRESS_FACTOR),
    gp1(MEMB_ADHE_RANGE),
    gp1(THRESH_SP),
    gp1(LJ_THRESH),
    //gp1(MEMB_ADJ_CONN_NUM),
    gp1(THRESH_DEAD),
    gpa(NUM_ITR,"LOOP"),
    gpa(Ca_avg_time,"Ca2P_avg_time"),
    gp1(ca2p_init),
    gp1(ex_inert_init),
    gp1(IP3_init),
    gp1(gj_init),
    gp1(ATP_init),
    gp1(ext_stim_init),
    gp1(FAC_MAP),
    gp1(SW_THRESH),
    gpa(T_TURNOVER,"TURNOVER"),
    gp1(NUM_SC_INIT),
    gp1(USE_TRI_MEMB),
    gp1(agki_max),
    gp1(fac),
    //gp1(agki_max_fix),
    gp1(div_max),
    //gp1(Ca_ITR),
    gp1(Kpp),
    gp1(dp),
    gp1(kbc),
    gp1(Hb),
    gp1(Cout),
    gp1(kg),
    gp1(gamma),
    gp1(mu0),
    gp1(mu1),
    gp1(kmu),
    gp1(para_b),
    gp1(para_bb),
    gp1(para_k1),
    gp1(k_flux),
    gp1(beta_zero),
    gp1(CA_OUT),
    //gp1(leak_from_storage),
    gp1(thgra),
    gp1(delta_th),
    gp1(thpri),
    gp1(kpa),
    //gp1(Kgra),
    gp1(delta_K),
    gp1(Kpri),
    gp1(delta_I),
    gp1(iage_kitei),
    gp1(para_k2),
    gp1(H0),
    gp1(wd),
    gp1(ca2p_du),
    gp1(STIM11),
    gp1(Kaa),
    gp1(Da),
    gp1(AIR_STIM),
    gp1(AREA_GRID_ORIGINAL),
    //gp1(AREA_GRID),
    //gp1(ANX),gp1(ANY),gp1(ANZ),
    gp1(N2),gp1(N3),
    gp1(eps_m),
    //gp1(P_MEMB),
    gp1(DER_DER_CONST),
    gp1(K_TOTAL),
    gp1(K_DESMOSOME_RATIO),
    //gp1(K_DESMOSOME),
    gp1(Kspring),
    gp1(Kspring_d),

    gp1(KBEND),
    gp1(delta_R_coef),
    //gp1(delta_R),
    gp1(para_ljp2),
    gp1(Kspring_division),
    gp1(stoch_div_time_ratio),
    gp1(S2),
    gp1(accel_div),
    gp1(eps_kb),
    gp1(alpha_b),
    gp1(S0),
    gp1(eps_ks),
    gp1(accel_diff),
    gp1(alpha_k),
    gp1(eps_kk),
    gp1(S1),
    gp1(ubar),
    gp1(delta_lipid),
    gp1(delta_sig_r1),
    gp1(lipid_rel),
    gp1(lipid),
        gp1(delta_L_coef),
    //gp1(delta_L),
    gp1(ADHE_CONST),
    gp1(DISA_conn_num_thresh),
    gp1(eps_L),
    gp1(unpair_dist_coef),
        gp1(kb),
        gp1(DUR_ALIVE),
        gp1(DUR_DEAD),
        gp1(DB),
        gp1(SYSTEM),
        gp1(CUT),
        gp1(NEW_BEND_POT),
        gp1(outputdir)
},
//set initial value
MEMB_NUM_X(100),
MEMB_NUM_Y(100),
MAX_CONNECT_CELL_NUM(400),
DT_Cell(0.01),
DT_Ca(0.01),
MALIG_NUM(0),
LX(50.0),LY(50.0),LZ(50.0),
NX(100),NY(100),NZ(100),
R_max(1.4),R_der(1.4),R_memb(1.0),
COMPRESS_FACTOR(4),
MEMB_ADHE_RANGE(0),
THRESH_SP(3.0),
LJ_THRESH(1.2),
THRESH_DEAD(22.0),
NUM_ITR(4* (unsigned int)1e6),
Ca_avg_time(10.0),
ca2p_init(0.122),
ex_inert_init(0.97),
IP3_init(0.0),
gj_init(0.99),
ATP_init(0.0),
ext_stim_init(0.0),
FAC_MAP(2.0),
SW_THRESH(20),
T_TURNOVER(6000.0),
NUM_SC_INIT(1),
USE_TRI_MEMB(false),
agki_max(6.0),
fac(1.0),
div_max(15),
Kpp(0.3),
dp(0.1),
kbc(0.4*1.2),
Hb(0.01),
Cout(1.0),
kg(0.1),
gamma(2.0),
mu0(0.567),
mu1(0.1),
kmu(0.05),
para_b(0.11),
para_bb(0.89),
para_k1(0.7),
k_flux(8.1),
beta_zero(0.02),
CA_OUT(1.0),
thgra(0.2),
delta_th(1.0),
thpri(1.0),
kpa(4.0),
delta_K(1.0),
Kpri(6.0),
delta_I(1.5),
iage_kitei(0.0),
para_k2(0.7),
H0(0.5),
wd(0.1),
ca2p_du(0.01),
STIM11(0.002),
Kaa(0.5),
Da(1.0),
AIR_STIM(0.1),
AREA_GRID_ORIGINAL(2.0),
N3(200),N2(400),
eps_m(0.01),
DER_DER_CONST(0.2),
K_TOTAL(3.0),
K_DESMOSOME_RATIO(0.01),
Kspring(25.0),
Kspring_d(5.0),
KBEND(0.5),
delta_R_coef(0.4),
para_ljp2(0.005),
Kspring_division(5.0),
stoch_div_time_ratio(0.25),
S2(0.1),
accel_div(1.0),
eps_kb(0.12),
alpha_b(5.0),
S0(0.1*0.2),
eps_ks(0.1*0.5),
accel_diff(1.0),
alpha_k(2.0),
eps_kk(0.1*0.5),
S1(0.1),
ubar(0.45),
delta_lipid(0.1),
delta_sig_r1(0.1),
lipid_rel(0.05*4.0),
lipid(0.6),
delta_L_coef(0.01),
ADHE_CONST(31.3),
DISA_conn_num_thresh(11),
eps_L(0.14),
unpair_dist_coef(0.9),
kb(0.025),
DUR_ALIVE(0.5),
DUR_DEAD(2.0),
DB(0.0009),
SYSTEM(0),
CUT(1000),
NEW_BEND_POT(true),
outputdir("output")
{
    calculate_param();
}

Params::Params(const std::string& path) :Params() {
    const_cast<Params*>(this)->load(path);
}
void Params::load(const std::string& path){
	auto kvmap = parse_paramtext(path);
    for (auto& pi : piset) {
        pi.set_from(kvmap);
    }
    calculate_param();
}
void Params::calculate_param() {
    agki_max_fix = fac*agki_max;
    Ca_ITR = (int)(Ca_avg_time / DT_Ca);
    leak_from_storage = CA_OUT*beta_zero;
    Kgra = kpa;
    AREA_GRID = AREA_GRID_ORIGINAL + 1e-7;
    ANX = (unsigned int)(LX / AREA_GRID_ORIGINAL + 0.5);
    ANY = (unsigned int)(LY / AREA_GRID_ORIGINAL + 0.5);
    ANZ = (unsigned int)(LZ / AREA_GRID_ORIGINAL);
    P_MEMB = 1.0 / COMPRESS_FACTOR;
    K_DESMOSOME = K_TOTAL*K_DESMOSOME_RATIO;
    delta_R = delta_R_coef*R_der;
    delta_L = delta_L_coef*R_max;

    MEMB_ADJ_CONN_NUM = USE_TRI_MEMB ? 6 : 4;
}

void Params::generate_paramfile(const char* out_path)const {
    std::ofstream of(out_path);
    if (!of) {
        throw std::runtime_error("Failed to create the file:"_s + out_path);
    }

    for (auto& pi : piset) {
        of << pi.name[0] << delim << pi.get_as_string()<<std::endl;
    }


}


