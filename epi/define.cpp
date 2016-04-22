#include <math.h>
#include "define.h"
//ファイルから後で読み込まれる
int EPI::C::NMEMB = -1;
int EPI::C::NDER = -1;
__m256d EPI::PState4d::all1 = _mm256_castsi256_pd(_mm256_set1_epi32(0xffffffff));
__m256d EPI::PState4d::all0 = _mm256_castsi256_pd(_mm256_set1_epi32(0x00000000));
VEC_INIT(EPI::C,zero);
VEC_INIT(EPI::C, half);
VEC_INIT(EPI::C, one);
VEC_INIT(EPI::C, neg_zero);
VEC_INIT(EPI::C,tanh_s_th);
VEC_INIT(EPI::C,tanh_l_th);
VEC_INIT(EPI::C,tanh_m_th);
VEC_INIT(EPI::C,tanh_dbl_p0);
VEC_INIT(EPI::C,tanh_dbl_p1);
VEC_INIT(EPI::C,tanh_dbl_p2);
VEC_INIT(EPI::C,tanh_dbl_q0);
VEC_INIT(EPI::C,tanh_dbl_q1);
VEC_INIT(EPI::C,tanh_dbl_q2);
VEC_INIT(EPI::C,tanh_dbl_q3);
VEC_INIT(EPI::C, dt_cell);
VEC_INIT(EPI::C, dt_ca);

VEC_INIT(EPI::C,Lx);
VEC_INIT(EPI::C,Ly);
VEC_INIT(EPI::C,Lz);
VEC_INIT(EPI::C, NX);
VEC_INIT(EPI::C, NY);
VEC_INIT(EPI::C, NZ);
VEC_INIT(EPI::C, dx);
VEC_INIT(EPI::C, dy);
VEC_INIT(EPI::C, dz);
VEC_INIT(EPI::C, dxSq);
VEC_INIT(EPI::C, dySq);
VEC_INIT(EPI::C, dzSq);

VEC_INIT(EPI::Ca2P,dA);
VEC_INIT(EPI::Ca2P,dP);
VEC_INIT(EPI::Ca2P,dc);
VEC_INIT(EPI::Ca2P,dB);
VEC_INIT(EPI::Ca2P,Kaa);
VEC_INIT(EPI::Ca2P,Kpp);
VEC_INIT(EPI::Ca2P,Kbb);
VEC_INIT(EPI::Ca2P,Kac);
VEC_INIT(EPI::Ca2P,Kf);
VEC_INIT(EPI::Ca2P,Kmu);
VEC_INIT(EPI::Ca2P,K1);
VEC_INIT(EPI::Ca2P,Kg);
VEC_INIT(EPI::Ca2P,Kbc);
VEC_INIT(EPI::Ca2P,mu0);
VEC_INIT(EPI::Ca2P,mu1);
VEC_INIT(EPI::Ca2P,alpha0);
VEC_INIT(EPI::Ca2P,sub1Alpha0);
VEC_INIT(EPI::Ca2P,gamma);
VEC_INIT(EPI::Ca2P,beta_zero);
VEC_INIT(EPI::Ca2P,CA_OUT);
VEC_INIT(EPI::Ca2P,beta);
VEC_INIT(EPI::Ca2P,Hb);
VEC_INIT(EPI::Ca2P,H0);
VEC_INIT(EPI::Ca2P,K2);
VEC_INIT(EPI::Ca2P,wd0);
VEC_INIT(EPI::Ca2P,eps_w0);
VEC_INIT(EPI::Ca2P,THRESH_DEAD);
VEC_INIT(EPI::Ca2P,DUR_DEAD);
VEC_INIT(EPI::Ca2P,DUR_ALIVE);
VEC_INIT(EPI::Ca2P,THRESH_SP);
VEC_INIT(EPI::Ca2P,Sk1);
VEC_INIT(EPI::Ca2P,Sk2);
VEC_INIT(EPI::Ca2P,Ss);
VEC_INIT(EPI::Ca2P,tau_g);
VEC_INIT(EPI::Ca2P,tau_s);
VEC_INIT(EPI::Ca2P,delta_tau);
VEC_INIT(EPI::Ca2P,delta_I);
VEC_INIT(EPI::Ca2P,delta_k);
VEC_INIT(EPI::Ca2P,kg);
VEC_INIT(EPI::Ca2P,ks);
VEC_INIT(EPI::Ca2P,_rSq);
void init() {
	for (int i = 0; i < EPI::C::NX; i++) {
		EPI::C::x_prev_arr[i]= ((i - 4 + EPI::C::NX) % EPI::C::NX) / 4;
		EPI::C::x_next_arr[i] = ((i + 4 + EPI::C::NX) % EPI::C::NX) / 4;
	}
	for (int i = 0; i < EPI::C::NY; i++) {
		EPI::C::y_prev_arr[i] = (i - 1 + EPI::C::NY) % EPI::C::NY;
		EPI::C::y_next_arr[i] = (i + 1 + EPI::C::NY) % EPI::C::NY;
	}
	////////////??
}

//for future use
__m256d _tanh_poly(const __m256d& v){
    __m256d sq = _mm256_mul_pd(v,v);
    __m256d __R = _mm256_div_pd(
                _mm256_mul_pd(
                    sq,
                    _mm256_add_pd(EPI::C::tanh_dbl_p0_4d,
                                  _mm256_mul_pd(sq,
                                                _mm256_add_pd(
                                                    EPI::C::tanh_dbl_p1_4d,
                                                    _mm256_mul_pd(sq,EPI::C::tanh_dbl_p2_4d)
                                                    )
                                                )
                                  )
                    ),
                _mm256_add_pd(EPI::C::tanh_dbl_q0_4d,
                              _mm256_mul_pd(sq,
                                            _mm256_add_pd(EPI::C::tanh_dbl_q1_4d,
                                                          _mm256_mul_pd(sq,
                                                                        _mm256_add_pd(EPI::C::tanh_dbl_q2_4d,
                                                                                      _mm256_mul_pd(sq,EPI::C::tanh_dbl_q3_4d)
                                                                                      )
                                                                        )
                                                          )
                                            )
                              )

                );
    //fmadd
    return _mm256_add_pd(v,_mm256_mul_pd(v,__R));
}
__m256d tanh_avx(const __m256d &x){
    //not optimized
    alignas(32) double stored[4];
    _mm256_store_pd(stored,x);
    stored[0]=tanh(stored[0]);
    stored[1]=tanh(stored[1]);
    stored[2]=tanh(stored[2]);
    stored[3]=tanh(stored[3]);
    return _mm256_load_pd(stored);
}

__m256d tanh_alt(const __m256d&x){
    //inefficient?
    return _mm256_mul_pd(EPI::C::half_4d,_mm256_add_pd(EPI::C::one_4d,tanh_avx(x)));
}

void VSet4d::operator+=(const VSet4d &a) {
	x = _mm256_add_pd(x, a.x);
	y = _mm256_add_pd(y, a.y);
	z = _mm256_add_pd(z, a.z);
}



namespace EPI{

	template<STATE, STATE...>
	__m256d PState4d::getORMask(STATE first, STATE... rest) {
		return _mm256_or_pd(smask[first], getORMask(rest...));
	}
	__m256d PState4d::getORMask() {
		return PState4d::all0;
	}

	template<STATE, STATE...>
	__m256d PState4d::getANDMask(STATE first, STATE... rest) {
		return _mm256_and_pd(smask[first], getANDMask(rest...));
	}
	__m256d PState4d::getANDMask() {
		return PState4d::all1;
	}
	void CellSet4d::get_lattice(VSet4d& out) const {
		__m256d raw_x_l = _mm256_floor_pd(_mm256_div_pd(pos.x, C::dx_4d));
		__m256d raw_y_l = _mm256_floor_pd(_mm256_div_pd(pos.y, C::dy_4d));
		__m256d raw_z_l = _mm256_floor_pd(_mm256_div_pd(pos.z, C::dz_4d));
		__m256d q_x = _mm256_floor_pd(_mm256_div_pd(raw_x_l, C::NX_4d));
		__m256d q_y = _mm256_floor_pd(_mm256_div_pd(raw_y_l, C::NY_4d));
		__m256d q_z = _mm256_floor_pd(_mm256_div_pd(raw_z_l, C::NZ_4d));
		out.x = _mm256_round_pd(_mm256_sub_pd(raw_x_l, _mm256_mul_pd(q_x, C::NX_4d)), _MM_FROUND_NINT);
		out.y = _mm256_round_pd(_mm256_sub_pd(raw_y_l, _mm256_mul_pd(q_y, C::NY_4d)), _MM_FROUND_NINT);
		out.z = _mm256_round_pd(_mm256_sub_pd(raw_z_l, _mm256_mul_pd(q_z, C::NZ_4d)), _MM_FROUND_NINT);
	}
__m256d Ca2P::G(const VSet4d& x_vec, const VSet4d& xi_vec, const __m256d& dcdt) {
    __m256d x_dist = _mm256_sub_pd(x_vec.x, xi_vec.x);
    __m256d y_dist = _mm256_sub_pd(x_vec.y, xi_vec.y);
    __m256d z_dist = _mm256_sub_pd(x_vec.z, xi_vec.z);
    x_dist = _mm256_mul_pd(x_dist, x_dist);
    y_dist = _mm256_mul_pd(y_dist, y_dist);
    z_dist = _mm256_mul_pd(z_dist, z_dist);

    //__m256 distSq = _mm256_add_ps(_mm256_add_ps(x_dist, y_dist), z_dist);
    /*
    __m256 kai=_mm256_cmp_ps(
        _mm256_add_ps(_mm256_add_ps(x_dist, y_dist), z_dist), //distSq
        _rSq8, _CMP_LE_OS);
    __m256 __R = _mm256_cmp_ps(dcdt, C::zero8, _CMP_GE_OS);
    */

    //X_0 && R(dc/dt) && K_ac
    return _mm256_and_pd(
                //mask
                _mm256_and_pd(
                    _mm256_cmp_pd(
                        _mm256_add_pd(_mm256_add_pd(x_dist, y_dist), z_dist),_rSq_4d, _CMP_LE_OS), //kai (|x-y|^2<=r^2)
                    _mm256_cmp_pd(dcdt, C::zero_4d, _CMP_GE_OS) //R(dc/dt) (dc/dt>=0)
                    ),
                Kac_4d
                );
}

__m256d Ca2P::Fc(const __m256d& P, const __m256d& c, const __m256d& h, const __m256d& B) {
    //分母揃えて一気に割った方が良い
    __m256d BSq = _mm256_mul_pd(B,B);

    __m256d den_Km_p = _mm256_add_pd(Kmu_4d,P);
    __m256d den_K1_c = _mm256_add_pd(K1_4d,c);
    __m256d den_Kg_c = _mm256_add_pd(Kg_4d,c);
    __m256d den_Hb_BSq = _mm256_add_pd(Hb_4d,BSq);

    __m256d num_mu1p = _mm256_mul_pd(mu1_4d,P);
    __m256d num_salphac = _mm256_mul_pd(sub1Alpha0_4d,c);
    __m256d num_KbcBSq = _mm256_mul_pd(Kbc_4d,BSq);
    __m256d num_gammac = _mm256_mul_pd(gamma_4d,c);

    //fmadd available if AVX2
    __m256d num_mu0Kmu_p = _mm256_add_pd(_mm256_mul_pd(mu0_4d,den_Km_p),num_mu1p);
    __m256d num_a0K1_c_salc = _mm256_add_pd(_mm256_mul_pd(alpha0_4d,den_K1_c),num_salphac);
    __m256d num_KbcBSqKg_c__gcHb_BSq = _mm256_sub_pd(_mm256_mul_pd(num_KbcBSq,den_Kg_c),_mm256_mul_pd(num_gammac,den_Hb_BSq));

    __m256d den_p1 = _mm256_mul_pd(den_Km_p,den_K1_c);
    __m256d den_p2 = _mm256_mul_pd(den_Kg_c,den_Hb_BSq);
    __m256d num_p1 = _mm256_mul_pd(_mm256_mul_pd(den_p2,_mm256_mul_pd(Kf_4d,h)),_mm256_mul_pd(num_mu0Kmu_p,num_a0K1_c_salc));
    __m256d num_p2 = _mm256_mul_pd(den_p1,num_KbcBSqKg_c__gcHb_BSq);
    return _mm256_add_pd(
                _mm256_div_pd(_mm256_add_pd(num_p1,num_p2),_mm256_mul_pd(den_p1,den_p2)),
                beta_4d);
    //add(sub):8 mul: 16
}
//inline?
__m256d Ca2P::FP(const __m256d& A,const __m256d& Ski){
    return _mm256_div_pd(_mm256_mul_pd(Kpa(Ski),A),_mm256_add_pd(H0_4d,A));
}
__m256d Ca2P::Fh(const __m256d &c, const __m256d & h){
    __m256d K2Sq = _mm256_mul_pd(K2_4d,K2_4d);
    return _mm256_sub_pd(_mm256_div_pd(K2Sq,_mm256_add_pd(c,K2Sq)),h);
}

__m256d Ca2P::Fw(const __m256d &wij, const __m256d &ci, const __m256d &cj){
    return _mm256_sub_pd(
                tanh_alt(
                    _mm256_div_pd(
                        _mm256_sub_pd(wd0_4d,_mm256_andnot_pd(neg_zero_4d,_mm256_sub_pd(ci,cj))), //calc'ing abs
                        eps_w0_4d
                        )
                    ),
                wij
                );
}
__m256d Ca2P::tau_h(const __m256d & Ski){
    return _mm256_add_pd(tau_g_4d,_mm256_mul_pd(_mm256_sub_pd(tau_s_4d,tau_g_4d),tanh_alt(_mm256_div_pd(_mm256_sub_pd(Ss_4d,Ski),delta_tau_4d))));
}
__m256d Ca2P::In(const __m256d &Ski){
    return tanh_alt(_mm256_div_pd(_mm256_sub_pd(Ski,Ss_4d),delta_I_4d));
}
__m256d Ca2P::Kpa(const __m256d &Ski){
    return _mm256_add_pd(kg_4d,_mm256_mul_pd(_mm256_sub_pd(ks_4d,kg_4d),tanh_alt(_mm256_div_pd(_mm256_sub_pd(Ss_4d,Ski),delta_k_4d))));
}

void Ca2P::refresh_Ca(const CUBE& calc_area,
	const _3DScalar4d& currentCa,
	const std::vector<CellSet4d>& all_cells,
	const std::vector<__m256d>& d_ci_dt,
	_3DScalar4d& nextCa) {
	//x is compressed by 4
	int prev_x_idx, next_x_idx, prev_y_idx, next_y_idx, prev_z_idx, next_z_idx;
	
	int _i,_i1,_i2,_i3,_i4;
	__m256d _sum;
	__m256d x_next_sub,x_prev_sub,y_next_sub,y_prev_sub,z_next_sub,z_prev_sub,sub_sum;
	__m256d x_tmp;
	__m256d dxSq;
	alignas(32) double x_tmp[8];
	VSet4d CaPos;
	//currently boundary conditon is not considered
	for (int i = calc_area.x.min; i <= calc_area.x.max/4; i++) {
		_i = i * 4;
		_i1 = _i*dx; _i2 = _i1 + dx; _i3 = _i2 + dx; _i4 = _i3 + dx;
		CaPos.x = _mm256_set_pd(_i1, _i2,_i3,_i4);
		//periodic boundary condition for x and y
		prev_x_idx = ((i - 4 + NX) % NX) / 4;
		next_x_idx = ((i + 4 + NX) % NX) / 4;
		for (int j = calc_area.y.min; j <= calc_area.y.max; j++) {
			CaPos.y = _mm256_set1_pd(i*dy);
			prev_y_idx = (j - 1 + NY) % NY;
			next_y_idx = (j + 1 + NY) % NY;
			for (int k = calc_area.z.min; k <= calc_area.z.max; k++) {
				CaPos.z = _mm256_set1_pd(i*dz);
				//neumann
				prev_z_idx = k == 0 ? 1 : k - 1;
				next_z_idx = k == NZ ? NZ - 1 : k + 1;
				for (int l = CELL_START_IDX; l < current_cell_num;l++){
					_sum = _mm256_add_pd(_sum, G(CaPos, all_cells[l].pos, d_ci_dt[l]));
				}
				//there is a more efficient way
				x_tmp = _mm256_set_pd(currentCa[i][j][k].m256d_f64[1], currentCa[i][j][k].m256d_f64[2], currentCa[i][j][k].m256d_f64[3], currentCa[next_x_idx][j][k].m256d_f64[0]);
				x_next_sub = _mm256_sub_pd(x_tmp, currentCa[i][j][k]);
				x_tmp = _mm256_set_pd(currentCa[prev_x_idx][j][k].m256d_f64[3], currentCa[i][j][k].m256d_f64[0], currentCa[i][j][k].m256d_f64[1], currentCa[i][j][k].m256d_f64[2]);
				x_prev_sub = _mm256_sub_pd(x_tmp, currentCa[i][j][k]);
				
				y_next_sub = _mm256_sub_pd(currentCa[i][next_y_idx][k], currentCa[i][j][k]);
				y_prev_sub = _mm256_sub_pd(currentCa[i][prev_y_idx][k], currentCa[i][j][k]);
				z_next_sub = _mm256_sub_pd(currentCa[i][j][next_z_idx], currentCa[i][j][k]);
				z_prev_sub = _mm256_sub_pd(currentCa[i][j][prev_z_idx], currentCa[i][j][k]);
				sub_sum = _mm256_add_pd(x_next_sub, _mm256_add_pd(x_prev_sub, _mm256_add_pd(y_next_sub, _mm256_add_pd(y_prev_sub, _mm256_add_pd(z_next_sub, z_prev_sub)))));
				nextCa[i][j][k] = _mm256_add_pd(
					currentCa[i][j][k],
					_mm256_mul_pd(dt_ca_4d,
						_mm256_add_pd(
							_mm256_mul_pd(dA_4d,_mm256_div_pd(sub_sum,dxSq_4d)),_mm256_sub_pd(_sum,_mm256_mul_pd(Kaa_4d, currentCa[i][j][k]))
							)
						)
					);

			}
		}
	}
}

void Ca2P::refresh_P_i(int calc_index_min, int calc_index_max,
	const _3DScalar4d& currentCa,
	const std::vector<CellSet4d>& all_current_cells,
	std::vector<CellSet4d>& refreshed_cells) {
	VSet4d lat;
	__m256d react_sum;
	__m256d neg_me_dead_val, react_alive_val, me_val, react_val;//dt unmultiplied
	for (int i = calc_index_min/4; i < calc_index_max/4; i++) {
		neg_me_dead_val = _mm256_and_pd(all_current_cells[i].state.smask[PState4d::DEAD],_mm256_mul_pd(Kpp_4d, all_current_cells[i].P));
		for (int j = 0; j < C::max_cell_num / 4; j++) {
			if (all_current_cells[i].react_flag_4d[j]) {
				react_alive_val = _mm256_add_pd(
					react_alive_val, 
					_mm256_and_pd(
						all_current_cells[j].state.smask[PState4d::ALIVE],
						_mm256_sub_pd(all_current_cells[j].P, all_current_cells[i].P)
						)
					);
			}
		}
		react_alive_val = _mm256_and_pd(all_current_cells[i].state.smask[PState4d::DEAD], _mm256_mul_pd(dP_4d, react_alive_val));
		all_current_cells[i].get_lattice(lat);
		//バラバラに計算しなきゃダメそう
		/*
		me_val = _mm256_and_pd(
			all_current_cells[i].state.getORMask(PState4d::ALIVE, PState4d::FIX, PState4d::MUSUME),FP()
		refreshed_cells[i].P = _mm256_add_pd()
		*/
		if(all_current_cells[i].state.DEAD)
		for (int j = 0; j < CellSet4d::max_reaction_cell_num; j++) {
			react_sum = _mm256_add_pd(react_sum,_mm256_and_pd(all_current_cells[i].react_mask[j],)
		}
	}
}
}
