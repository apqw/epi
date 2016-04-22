#include <math.h>
#include "define.h"
//?t?@?C??????????????????
int EPI::C::NMEMB = -1;
int EPI::C::NDER = -1;
int EPI::C::CELL_START_IDX = NMEMB + NDER; //tmp def
__m256d EPI::PState4d::all1 = _mm256_castsi256_pd(_mm256_set1_epi32(0xffffffff));
__m256d EPI::PState4d::all0 = _mm256_castsi256_pd(_mm256_set1_epi32(0x00000000));
int EPI::Field_Data::current_cell_num = 0;
std::vector<EPI::CellSet4d> EPI::Field_Data::cells(C::max_cell_num/4);
std::vector<EPI::CellSet4d> EPI::Field_Data::next_cells(C::max_cell_num/4);
//int EPI::current_cell_num = 0;

VEC_INIT(EPI::C, zero);
VEC_INIT(EPI::C, half);
VEC_INIT(EPI::C, one);
VEC_INIT(EPI::C, neg_zero);
VEC_INIT(EPI::C, tanh_s_th);
VEC_INIT(EPI::C, tanh_l_th);
VEC_INIT(EPI::C, tanh_m_th);
VEC_INIT(EPI::C, tanh_dbl_p0);
VEC_INIT(EPI::C, tanh_dbl_p1);
VEC_INIT(EPI::C, tanh_dbl_p2);
VEC_INIT(EPI::C, tanh_dbl_q0);
VEC_INIT(EPI::C, tanh_dbl_q1);
VEC_INIT(EPI::C, tanh_dbl_q2);
VEC_INIT(EPI::C, tanh_dbl_q3);
VEC_INIT(EPI::C, dt_cell);
VEC_INIT(EPI::C, dt_ca);

VEC_INIT(EPI::C, Lx);
VEC_INIT(EPI::C, Ly);
VEC_INIT(EPI::C, Lz);
VEC_INIT(EPI::C, NX);
VEC_INIT(EPI::C, NY);
VEC_INIT(EPI::C, NZ);
VEC_INIT(EPI::C, dx);
VEC_INIT(EPI::C, dy);
VEC_INIT(EPI::C, dz);
VEC_INIT(EPI::C, dxSq);
VEC_INIT(EPI::C, dySq);
VEC_INIT(EPI::C, dzSq);

VEC_INIT(EPI::Ca2P, dA);
VEC_INIT(EPI::Ca2P, dP);
VEC_INIT(EPI::Ca2P, dc);
VEC_INIT(EPI::Ca2P, dB);
VEC_INIT(EPI::Ca2P, Kaa);
VEC_INIT(EPI::Ca2P, Kpp);
VEC_INIT(EPI::Ca2P, Kbb);
VEC_INIT(EPI::Ca2P, Kac);
VEC_INIT(EPI::Ca2P, Kf);
VEC_INIT(EPI::Ca2P, Kmu);
VEC_INIT(EPI::Ca2P, K1);
VEC_INIT(EPI::Ca2P, Kg);
VEC_INIT(EPI::Ca2P, Kbc);
VEC_INIT(EPI::Ca2P, mu0);
VEC_INIT(EPI::Ca2P, mu1);
VEC_INIT(EPI::Ca2P, alpha0);
VEC_INIT(EPI::Ca2P, sub1Alpha0);
VEC_INIT(EPI::Ca2P, gamma);
VEC_INIT(EPI::Ca2P, beta_zero);
VEC_INIT(EPI::Ca2P, CA_OUT);
VEC_INIT(EPI::Ca2P, beta);
VEC_INIT(EPI::Ca2P, Hb);
VEC_INIT(EPI::Ca2P, H0);
VEC_INIT(EPI::Ca2P, K2);
VEC_INIT(EPI::Ca2P, wd0);
VEC_INIT(EPI::Ca2P, eps_w0);
VEC_INIT(EPI::Ca2P, THRESH_DEAD);
VEC_INIT(EPI::Ca2P, DUR_DEAD);
VEC_INIT(EPI::Ca2P, DUR_ALIVE);
VEC_INIT(EPI::Ca2P, THRESH_SP);
VEC_INIT(EPI::Ca2P, Sk1);
VEC_INIT(EPI::Ca2P, Sk2);
VEC_INIT(EPI::Ca2P, Ss);
VEC_INIT(EPI::Ca2P, tau_g);
VEC_INIT(EPI::Ca2P, tau_s);
VEC_INIT(EPI::Ca2P, delta_tau);
VEC_INIT(EPI::Ca2P, delta_I);
VEC_INIT(EPI::Ca2P, delta_k);
VEC_INIT(EPI::Ca2P, kg);
VEC_INIT(EPI::Ca2P, ks);
VEC_INIT(EPI::Ca2P, _rSq);
VEC_INIT(EPI::Ca2P, iage_kitei);
VEC_INIT(EPI::Ca2P, Cout);
void init() {
	/*
	for (int i = 0; i < EPI::C::NX; i++) {
		EPI::C::x_prev_arr[i]= ((i - 4 + EPI::C::NX) % EPI::C::NX) / 4;
		EPI::C::x_next_arr[i] = ((i + 4 + EPI::C::NX) % EPI::C::NX) / 4;
	}
	for (int i = 0; i < EPI::C::NY; i++) {
		EPI::C::y_prev_arr[i] = (i - 1 + EPI::C::NY) % EPI::C::NY;
		EPI::C::y_next_arr[i] = (i + 1 + EPI::C::NY) % EPI::C::NY;
	}
	////////////??
	/// */
}

//for future use
__m256d _tanh_poly(const __m256d& v) {
	__m256d sq = _mm256_mul_pd(v, v);
	__m256d __R = _mm256_div_pd(
		_mm256_mul_pd(
			sq,
			_mm256_add_pd(EPI::C::tanh_dbl_p0_4d,
				_mm256_mul_pd(sq,
					_mm256_add_pd(
						EPI::C::tanh_dbl_p1_4d,
						_mm256_mul_pd(sq, EPI::C::tanh_dbl_p2_4d)
						)
					)
				)
			),
		_mm256_add_pd(EPI::C::tanh_dbl_q0_4d,
			_mm256_mul_pd(sq,
				_mm256_add_pd(EPI::C::tanh_dbl_q1_4d,
					_mm256_mul_pd(sq,
						_mm256_add_pd(EPI::C::tanh_dbl_q2_4d,
							_mm256_mul_pd(sq, EPI::C::tanh_dbl_q3_4d)
							)
						)
					)
				)
			)

		);
	//fmadd
	return _mm256_add_pd(v, _mm256_mul_pd(v, __R));
}
__m256d tanh_avx(const __m256d &x) {
	//not optimized
	alignas(32) double stored[4];
	_mm256_store_pd(stored, x);
	stored[0] = tanh(stored[0]);
	stored[1] = tanh(stored[1]);
	stored[2] = tanh(stored[2]);
	stored[3] = tanh(stored[3]);
	return _mm256_load_pd(stored);
}

__m256d tanh_alt(const __m256d&x) {
	//inefficient?
	return _mm256_mul_pd(EPI::C::half_4d, _mm256_add_pd(EPI::C::one_4d, tanh_avx(x)));
}
__m256d m256dintmod(const __m256d &num, const __m256d& den) {
	return _mm256_round_pd(
		_mm256_sub_pd(num,
			_mm256_mul_pd(
				_mm256_floor_pd(
					_mm256_div_pd(num, den)
					), den
				)
			),
		_MM_FROUND_NINT);
	//Round(num-Floor(num/den)*den)
}

void VSet4d::operator+=(const VSet4d &a) {
	x = _mm256_add_pd(x, a.x);
	y = _mm256_add_pd(y, a.y);
	z = _mm256_add_pd(z, a.z);
}

__m256d calc_avg8(const VSet4d& lat, const _3DScalar4d& _3DVal_4d) {
	//平均をとる番号
	std::vector<VSet4d> avg_lat(8);
	alignas(16) int avg_index_y[4], avg_index_z[4],avg_vec_index_x[4], avg_invec_index_x[4];
	alignas(32) double avg_vec_store[4],avg_val[4];
	__m256d avg_vec_idx_d = { 0 }, avg = { 0 };
	avg_lat[0].x = lat.x;
	avg_lat[0].y = lat.y;
	avg_lat[0].z = lat.z;

	//boundary is not considered
	avg_lat[1].x = _mm256_add_pd(_mm256_set1_pd(1.0), lat.x);
	avg_lat[1].y = lat.y;
	avg_lat[1].z = lat.z;

	avg_lat[2].x = lat.x;
	avg_lat[2].y = _mm256_add_pd(_mm256_set1_pd(1.0), lat.y);
	avg_lat[2].z = lat.z;

	avg_lat[3].x = lat.x;
	avg_lat[3].y = lat.y;
	avg_lat[3].z = _mm256_add_pd(_mm256_set1_pd(1.0), lat.z);

	avg_lat[4].x = lat.x;
	avg_lat[4].y = avg_lat[2].y; //reuse that is already added
	avg_lat[4].z = avg_lat[3].z;

	avg_lat[5].x = avg_lat[1].x;
	avg_lat[5].y = lat.y;
	avg_lat[5].z = avg_lat[3].z;

	avg_lat[6].x = avg_lat[1].x;
	avg_lat[6].y = avg_lat[2].y;
	avg_lat[6].z = lat.z;

	avg_lat[7].x = avg_lat[1].x;
	avg_lat[7].y = avg_lat[2].y;
	avg_lat[7].z = avg_lat[3].z;

	for (int k = 0; k < 8; k++) {

		_mm_store_si128((__m128i*)avg_index_y, _mm256_cvtpd_epi32(avg_lat[k].y));
		_mm_store_si128((__m128i*)avg_index_z, _mm256_cvtpd_epi32(avg_lat[k].z));
		avg_vec_idx_d = _mm256_round_pd(_mm256_div_pd(avg_lat[k].x, _mm256_set1_pd(4.0)), _MM_FROUND_NINT); //index as double
		_mm_store_si128((__m128i*)avg_vec_index_x, _mm256_cvtpd_epi32(avg_vec_idx_d)); //index of 4-pair data
		_mm_store_si128(
			(__m128i*)avg_invec_index_x,
			_mm256_cvtpd_epi32(
				_mm256_sub_pd(avg_lat[k].x,
					_mm256_mul_pd(avg_vec_idx_d, _mm256_set1_pd(4.0))
					)
				)
			); //index of element in 4-pairs [0,1,2,3]

			   //仕方ないので４つバラバラにして取得
		for (int cell_count = 0; cell_count < 4; cell_count++) {
			_mm256_store_pd(avg_vec_store, _3DVal_4d[avg_vec_index_x[cell_count]][avg_index_y[cell_count]][avg_index_z[cell_count]]);
			avg_val[cell_count] = avg_vec_store[avg_invec_index_x[cell_count]];
		}
		avg = _mm256_add_pd(avg, _mm256_load_pd(avg_val));
	}
	return _mm256_div_pd(avg, _mm256_set1_pd(8.0)); //4箇所それぞれの平均が入ってる

	/* original code
	
	avg_lat[0].x = lat.x;
	avg_lat[0].y = lat.y;
	avg_lat[0].z = lat.z;

	//boundary is not considered
	avg_lat[1].x = _mm256_add_pd(_mm256_set1_pd(1.0), lat.x);
	avg_lat[1].y = lat.y;
	avg_lat[1].z = lat.z;

	avg_lat[2].x = lat.x;
	avg_lat[2].y = _mm256_add_pd(_mm256_set1_pd(1.0), lat.y);
	avg_lat[2].z = lat.z;

	avg_lat[3].x = lat.x;
	avg_lat[3].y = lat.y;
	avg_lat[3].z = _mm256_add_pd(_mm256_set1_pd(1.0), lat.z);

	avg_lat[4].x = lat.x;
	avg_lat[4].y = avg_lat[2].y; //reuse that is already added
	avg_lat[4].z = avg_lat[3].z;

	avg_lat[5].x = avg_lat[1].x;
	avg_lat[5].y = lat.y;
	avg_lat[5].z = avg_lat[3].z;

	avg_lat[6].x = avg_lat[1].x;
	avg_lat[6].y = avg_lat[2].y;
	avg_lat[6].z = lat.z;

	avg_lat[7].x = avg_lat[1].x;
	avg_lat[7].y = avg_lat[2].y;
	avg_lat[7].z = avg_lat[3].z;

	for (int k = 0; k < 8; k++) {

	_mm_store_si128((__m128i*)avg_index_y, _mm256_cvtpd_epi32(avg_lat[k].y));
	_mm_store_si128((__m128i*)avg_index_z, _mm256_cvtpd_epi32(avg_lat[k].z));
	avg_vec_idx_d = _mm256_round_pd(_mm256_div_pd(avg_lat[k].x, _mm256_set1_pd(4.0)), _MM_FROUND_NINT); //index as double
	_mm_store_si128((__m128i*)avg_vec_index_x, _mm256_cvtpd_epi32(avg_vec_idx_d)); //index of 4-pair data
	_mm_store_si128(
	(__m128i*)avg_invec_index_x,
	_mm256_cvtpd_epi32(
	_mm256_sub_pd(avg_lat[k].x,
	_mm256_mul_pd(avg_vec_idx_d, _mm256_set1_pd(4.0))
	)
	)
	); //index of element in 4-pairs [0,1,2,3]

	//仕方ないので４つバラバラにして取得
	for (int cell_count = 0; cell_count < 4; cell_count++) {
	_mm256_store_pd(avg_vec_store, currentCa[avg_vec_index_x[cell_count]][avg_index_y[cell_count]][avg_index_z[cell_count]]);
	avg_Ca_value[cell_count] = avg_vec_store[avg_invec_index_x[cell_count]];
	}
	avg = _mm256_add_pd(avg, _mm256_load_pd(avg_Ca_value));
	}
	avg = _mm256_div_pd(avg, _mm256_set1_pd(8.0)); //4箇所それぞれの平均が入ってる

	
	
	*/
}


namespace EPI {

	template<typename T, typename... U>
	__m256d PState4d::getORMask(T first, U... rest) const {
		return _mm256_or_pd(smask[first], getORMask(rest...));
	}
	__m256d PState4d::getORMask() const {
		return PState4d::all0;
	}

	template<typename T, typename... U>
	__m256d PState4d::getANDMask(T first, U... rest) const {
		return _mm256_and_pd(smask[first], getANDMask(rest...));
	}
	__m256d PState4d::getANDMask() const {
		return PState4d::all1;
	}


	void CellSet4d::get_lattice(VSet4d& out) const {
		__m256d raw_x_l = _mm256_floor_pd(_mm256_div_pd(pos.x, C::dx_4d));
		__m256d raw_y_l = _mm256_floor_pd(_mm256_div_pd(pos.y, C::dy_4d));
		__m256d raw_z_l = _mm256_floor_pd(_mm256_div_pd(pos.z, C::dz_4d));
		out.x = m256dintmod(raw_x_l, C::NX_4d);
		out.y = m256dintmod(raw_y_l, C::NY_4d);
		out.z = m256dintmod(raw_z_l, C::NZ_4d);
	}
	template<typename T, typename... U>
	bool CellSet4d::hasState(T s, U... rest) const {
		return (_mm256_testz_pd(this->state.smask[s], this->state.smask[s]) == 0 || hasState(rest...));
	}

	bool CellSet4d::hasState()const {
		return false;
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
					_mm256_add_pd(_mm256_add_pd(x_dist, y_dist), z_dist), _rSq_4d, _CMP_LE_OS), //kai (|x-y|^2<=r^2)
				_mm256_cmp_pd(dcdt, C::zero_4d, _CMP_GE_OS) //R(dc/dt) (dc/dt>=0)
				),
			Kac_4d
			);
	}

	__m256d Ca2P::Fc(const CellSet4d& cset, const __m256d& B) {
		//??????????C??????????????
		__m256d BSq = _mm256_mul_pd(B, B);

		__m256d den_Km_p = _mm256_add_pd(Kmu_4d, cset.P);
		__m256d den_K1_c = _mm256_add_pd(K1_4d, cset.c);
		__m256d den_Kg_c = _mm256_add_pd(Kg_4d, cset.c);
		__m256d den_Hb_BSq = _mm256_add_pd(Hb_4d, BSq);

		__m256d num_mu1p = _mm256_mul_pd(mu1_4d, cset.P);
		__m256d num_salphac = _mm256_mul_pd(sub1Alpha0_4d, cset.c);
		__m256d num_KbcBSq = _mm256_mul_pd(_mm256_mul_pd(Kbc_4d,Cout_4d), BSq); //why cout?
		__m256d num_gammac = _mm256_mul_pd(gamma_4d, cset.c);

		//fmadd available if AVX2
		__m256d num_mu0Kmu_p = _mm256_add_pd(_mm256_mul_pd(mu0_4d, den_Km_p), num_mu1p);
		__m256d num_a0K1_c_salc = _mm256_add_pd(_mm256_mul_pd(alpha0_4d, den_K1_c), num_salphac);
		__m256d num_KbcBSqKg_c__gcHb_BSq = _mm256_sub_pd(_mm256_mul_pd(num_KbcBSq, den_Kg_c), _mm256_mul_pd(num_gammac, den_Hb_BSq));

		__m256d den_p1 = _mm256_mul_pd(den_Km_p, den_K1_c);
		__m256d den_p2 = _mm256_mul_pd(den_Kg_c, den_Hb_BSq);
		__m256d num_p1 = _mm256_mul_pd(_mm256_mul_pd(den_p2, _mm256_mul_pd(Kf_4d, cset.h)), _mm256_mul_pd(num_mu0Kmu_p, num_a0K1_c_salc));
		__m256d num_p2 = _mm256_mul_pd(den_p1, num_KbcBSqKg_c__gcHb_BSq);
		return _mm256_add_pd(
			_mm256_div_pd(_mm256_add_pd(num_p1, num_p2), _mm256_mul_pd(den_p1, den_p2)),
			beta_4d);
		//add(sub):8 mul: 17
	}
	//inline?
	__m256d Ca2P::FP(const CellSet4d& cset, const __m256d& A) {
		return _mm256_div_pd(_mm256_mul_pd(Kpa(cset), A), _mm256_add_pd(H0_4d, A));
	}
	__m256d Ca2P::Fh(const CellSet4d& cset) {
		__m256d K2Sq = _mm256_mul_pd(K2_4d, K2_4d);
		return _mm256_sub_pd(_mm256_div_pd(K2Sq, _mm256_add_pd(cset.c, K2Sq)), cset.h);
	}

	__m256d Ca2P::Fw(const __m256d &wij, const __m256d &ci, const __m256d &cj) {
		return _mm256_sub_pd(
			tanh_alt(
				_mm256_div_pd(
					_mm256_sub_pd(wd0_4d, _mm256_andnot_pd(neg_zero_4d, _mm256_sub_pd(ci, cj))), //calc'ing abs
					eps_w0_4d
					)
				),
			wij
			);
	}
	__m256d Ca2P::tau_h(const CellSet4d& cset) {

		__m256d mask = cset.state.getORMask(FIX, MUSUME);
		__m256d age = _mm256_or_pd(_mm256_and_pd(mask, cset.ageb), _mm256_andnot_pd(mask, cset.agek));
		__m256d alive_val = _mm256_add_pd(tau_g_4d, _mm256_mul_pd(_mm256_sub_pd(tau_s_4d, tau_g_4d), tanh_alt(_mm256_div_pd(_mm256_sub_pd(Ss_4d, age), delta_tau_4d))));
		//__m256d fixmusume_val = ks_4d;
		__m256d not_alive_but_fix_or_musume = _mm256_andnot_pd(cset.state.smask[ALIVE], mask);
		__m256d not_alive_nor_fix_nor_musume = _mm256_andnot_pd(cset.state.smask[ALIVE], _mm256_andnot_pd(mask, PState4d::all1));

		return _mm256_or_pd(

			_mm256_and_pd(cset.state.smask[ALIVE], alive_val),

			_mm256_or_pd(
				_mm256_and_pd(not_alive_but_fix_or_musume, tau_s_4d),
				_mm256_and_pd(not_alive_nor_fix_nor_musume, zero_4d)
				)
			);
		//return _mm256_add_pd(tau_g_4d, _mm256_mul_pd(_mm256_sub_pd(tau_s_4d, tau_g_4d), tanh_alt(_mm256_div_pd(_mm256_sub_pd(Ss_4d, Ski), delta_tau_4d))));
	}
	__m256d Ca2P::In(const CellSet4d& cset) {
		__m256d mask = cset.state.getORMask(FIX, MUSUME);

		//use agek?
		return _mm256_or_pd(_mm256_and_pd(mask, iage_kitei_4d), _mm256_andnot_pd(mask, tanh_alt(_mm256_div_pd(_mm256_sub_pd(cset.agek, Ss_4d), delta_I_4d))));
		//return tanh_alt(_mm256_div_pd(_mm256_sub_pd(Ski,Ss_4d),delta_I_4d));
	}
	__m256d Ca2P::Kpa(const CellSet4d& cset) {

		__m256d mask = cset.state.getORMask(FIX, MUSUME);
		__m256d age = _mm256_or_pd(_mm256_and_pd(mask, cset.ageb), _mm256_andnot_pd(mask, cset.agek));
		__m256d alive_val = _mm256_add_pd(kg_4d, _mm256_mul_pd(_mm256_sub_pd(ks_4d, kg_4d), tanh_alt(_mm256_div_pd(_mm256_sub_pd(Ss_4d, age), delta_k_4d))));
		//__m256d fixmusume_val = ks_4d;
		__m256d not_alive_but_fix_or_musume = _mm256_andnot_pd(cset.state.smask[ALIVE], mask);
		__m256d not_alive_nor_fix_nor_musume = _mm256_andnot_pd(cset.state.smask[ALIVE], _mm256_andnot_pd(mask, PState4d::all1));

		return _mm256_or_pd(

			_mm256_and_pd(cset.state.smask[ALIVE], alive_val),

			_mm256_or_pd(
				_mm256_and_pd(not_alive_but_fix_or_musume, ks_4d),
				_mm256_and_pd(not_alive_nor_fix_nor_musume, zero_4d) //intentionally cause zero div -> error
				)
			);


		//return _mm256_add_pd(kg_4d,_mm256_mul_pd(_mm256_sub_pd(ks_4d,kg_4d),tanh_alt(_mm256_div_pd(_mm256_sub_pd(Ss_4d,age),delta_k_4d))));
	}

	void Ca2P::refresh_Ca(const CUBE& calc_area,
		const _3DScalar4d& currentCa,
		const std::vector<CellSet4d>& all_cells,
		_3DScalar4d& nextCa) {
		//x is compressed by 4
		VEC_X_RANGE_VALIDATION(calc_area.x.min); VEC_X_RANGE_VALIDATION(calc_area.x.max);

		int prev_x_idx, next_x_idx, prev_y_idx, next_y_idx, prev_z_idx, next_z_idx;

		double _i, _i1, _i2, _i3, _i4;
		__m256d _sum = { 0 };
		__m256d x_next_sub, x_prev_sub, y_next_sub, y_prev_sub, z_next_sub, z_prev_sub, sub_sum;
		__m256d x_tmp;
		alignas(32) double x_store_prev[4], x_store_next[4], x_store_current[4];
		VSet4d CaPos;
		//currently boundary conditon is not considered
		for (int i = calc_area.x.min / 4; i <= calc_area.x.max / 4; i++) {
			_i = i * 4;
			_i1 = _i*dx; _i2 = _i1 + dx; _i3 = _i2 + dx; _i4 = _i3 + dx;
			CaPos.x = _mm256_set_pd(_i1, _i2, _i3, _i4);
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
					for (int l = CELL_START_IDX; l < Field_Data::current_cell_num; l++) {
						_sum = _mm256_add_pd(_sum, G(CaPos, all_cells[l].pos, all_cells[l].diff_c));
					}
					_mm256_store_pd(x_store_prev, currentCa[prev_x_idx][j][k]);
					_mm256_store_pd(x_store_next, currentCa[next_x_idx][j][k]);
					_mm256_store_pd(x_store_current, currentCa[i][j][k]);
					//there is a more efficient way
					x_tmp = _mm256_set_pd(x_store_current[1], x_store_current[2], x_store_current[3], x_store_next[0]);
					x_next_sub = _mm256_sub_pd(x_tmp, currentCa[i][j][k]);
					x_tmp = _mm256_set_pd(x_store_prev[3], x_store_current[0], x_store_current[1], x_store_current[2]);
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
								_mm256_mul_pd(dA_4d, _mm256_div_pd(sub_sum, dxSq_4d)), _mm256_sub_pd(_sum, _mm256_mul_pd(Kaa_4d, currentCa[i][j][k]))
								)
							)
						);

				}
			}
		}
	}

	//indices must be 4n-aligned
	void Ca2P::refresh_P_i(int calc_index_min, int calc_index_max, \
		const _3DScalar4d& currentCa, \
		const std::vector<CellSet4d>& all_current_cells, \
		std::vector<CellSet4d>& refreshed_cells) {

		VEC_X_RANGE_VALIDATION(calc_index_min);
		VEC_X_RANGE_VALIDATION(calc_index_max);

		VSet4d lat; int lat_x_val = 0;
		/*
		std::vector<VSet4d> avg_lat(8);
		alignas(32) double avg_vec_store[4] = { 0 }, avg_Ca_value[4] = { 0 };
		__m256d avg = { 0 };
		*/
		alignas(16) int avg_index_x[4] = { 0 }, avg_index_y[4] = { 0 }, avg_index_z[4] = { 0 }, avg_vec_index_x[4] = { 0 }, avg_invec_index_x[4] = { 0 };
		__m256d avg_vec_idx_d = { 0 };
		__m256d react_sum = { 0 };
		__m256d neg_me_dead_val = { 0 }, react_alive_val = { 0 }, me_val = { 0 }, react_val = { 0 };//dt unmultiplied
		__m256d first_calc_val = { 0 }, second_calc_val = { 0 };
		__m256d react_mask = { 0 };
		__m256d In_val = { 0 };
		__m256d main_cond_mask = { 0 };
		for (int i = calc_index_min / 4; i < calc_index_max / 4; i++) {
			neg_me_dead_val = { 0 };
			react_alive_val = { 0 };
			me_val = { 0 };
			react_val = { 0 };

			//first calc
			if (all_current_cells[i].hasState(DEAD)) {
				//Kpp*Pi
				neg_me_dead_val = _mm256_mul_pd(Kpp_4d, all_current_cells[i].P);
				for (int j = 0; j < C::max_cell_num / 4; j++) {
					if (all_current_cells[j].hasState(ALIVE)) { //これでは全ての細胞に対して計算する
						//やっぱり重いなら変える
						//Pj-Pi
						react_alive_val = _mm256_add_pd(
							react_alive_val,
							_mm256_and_pd(
								all_current_cells[j].state.smask[ALIVE],
								_mm256_sub_pd(all_current_cells[j].P, all_current_cells[i].P)
								)
							);
					}
				}
				//dp*Sum(Pj-Pi)
				react_alive_val = _mm256_mul_pd(dP_4d, react_alive_val);
			}
			first_calc_val = _mm256_and_pd(all_current_cells[i].state.smask[DEAD], _mm256_sub_pd(react_alive_val, neg_me_dead_val));
			//先にDEADマスクをかけずにneg_me_dead_val+react_alive_valを計算してからかけたほうが良い?
			//->変わらない,けど見やすいかも

			main_cond_mask = all_current_cells[i].state.getORMask(ALIVE, FIX, MUSUME);
			if (_mm256_testz_pd(main_cond_mask, main_cond_mask) == 0) {
				//座標を格子点の番号に変換
				all_current_cells[i].get_lattice(lat);

				//avg = calc_avg8(lat, currentCa);

				//FP(A_avg)-Kpp*Pi
				me_val = _mm256_sub_pd(FP(all_current_cells[i], calc_avg8(lat, currentCa)), _mm256_mul_pd(Kpp_4d, all_current_cells[i].P));

				//In(Ski)
				In_val = In(all_current_cells[i]);


				for (int j = 0; j < C::max_cell_num / 4; j++) {
					//とりあえず全部

					react_mask = all_current_cells[j].state.getORMask(ALIVE, DEAD, FIX, MUSUME);

					if (_mm256_testz_pd(react_mask, react_mask) == 0) { //対象がある hasState(ALIVE,DEAD,FIX,MUSUME)も可だがマスクを使うので

						//ここを他の計算と一緒にしたら早くなるけどとりあえず放置..読みにくい



						//In*w_ij*(Pj-Pi)  Pj-Piさっきやったから使えそう


						react_val = _mm256_add_pd(react_val, _mm256_and_pd(react_mask,
							_mm256_mul_pd(
								In_val,
								_mm256_mul_pd(all_current_cells[i].w[j],
									_mm256_sub_pd(all_current_cells[j].P, all_current_cells[i].P)
									)
								)
							)
							);

					}
				}
				//dp*Sum(In*w_ij*(Pj-Pi))
				react_val = _mm256_mul_pd(dP_4d, react_val);

			}

			second_calc_val = _mm256_and_pd(main_cond_mask, _mm256_add_pd(me_val, react_val));

			refreshed_cells[i].P = _mm256_add_pd(all_current_cells[i].P,
				_mm256_mul_pd(C::dt_ca_4d, _mm256_add_pd(first_calc_val, second_calc_val))
				);
		}
	}

	void Ca2P::refresh_c_i(int calc_index_min, int calc_index_max,
		const _3DScalar4d& currentB,
		const std::vector<CellSet4d>& all_current_cells,
		std::vector<__m256d>& diff_c_out,
		std::vector<CellSet4d>& refreshed_cells) {
		VEC_X_RANGE_VALIDATION(calc_index_min);
		VEC_X_RANGE_VALIDATION(calc_index_max);
		__m256d main_cond_mask = { 0 };
		__m256d react_mask = { 0 };
		__m256d In_val = { 0 }, react_val = { 0 };
		VSet4d lat;
		std::fill(diff_c_out.begin(), diff_c_out.end(), C::zero_4d);
		for (int i = calc_index_min / 4; i < calc_index_max / 4; i++) {
			react_val = { 0 };
			main_cond_mask = all_current_cells[i].state.getORMask(ALIVE, FIX, MUSUME);
			if (_mm256_testz_pd(main_cond_mask, main_cond_mask) == 0) {
				all_current_cells[i].get_lattice(lat);
				//default diffu=0
				diff_c_out[i] = Fc(all_current_cells[i], calc_avg8(lat, currentB)); //mask after!!
				In_val = In(all_current_cells[i]);
				for (int j = 0; j < C::max_cell_num / 4; j++) {
					react_mask = all_current_cells[j].state.getORMask(ALIVE, DEAD, FIX, MUSUME);
					if (_mm256_testz_pd(react_mask, react_mask) == 0) {
						react_val = _mm256_add_pd(react_val, _mm256_and_pd(react_mask,
							_mm256_mul_pd(
								In_val,
								_mm256_mul_pd(all_current_cells[i].w[j],
									_mm256_sub_pd(all_current_cells[j].c, all_current_cells[i].c)
									)
								)
							)
							);
					}
				}

				//mask over all value,
				diff_c_out[i] = _mm256_and_pd(main_cond_mask,_mm256_add_pd(diff_c_out[i], _mm256_add_pd(dc_4d, react_val)));
				
				
				
				
			}
			//0 when the cond not satisfied
			//see -> u_ave[j]=u[j]=0.0;
			refreshed_cells[i].c = _mm256_and_pd(main_cond_mask, _mm256_add_pd(all_current_cells[i].c, _mm256_mul_pd(C::dt_ca_4d, diff_c_out[i])));
		}
	}

	void Ca2P::refresh_h_i(int calc_index_min, int calc_index_max,
		const std::vector<CellSet4d>& all_current_cells,
		std::vector<CellSet4d>& refreshed_cells) {
		VEC_X_RANGE_VALIDATION(calc_index_min);
		VEC_X_RANGE_VALIDATION(calc_index_max);

		__m256d main_cond_mask = { 0 };
		__m256d me_val = { 0 };
		for (int i = calc_index_min / 4; i < calc_index_max / 4; i++) {
			me_val = { 0 };
			main_cond_mask = all_current_cells[i].state.getORMask(ALIVE, FIX, MUSUME);
			if (_mm256_testz_pd(main_cond_mask, main_cond_mask) == 0) {

				me_val = _mm256_add_pd(me_val,
					_mm256_and_pd(main_cond_mask,
						_mm256_div_pd(Fh(all_current_cells[i]),tau_h(all_current_cells[i]))
						)
					);
			}
			refreshed_cells[i].h = _mm256_add_pd(
				all_current_cells[i].h, _mm256_and_pd(main_cond_mask, _mm256_and_pd(dt_ca_4d, me_val))
				);

		}
	}

	void Ca2P::refresh_w_i_j(int calc_index_min, int calc_index_max,
		const std::vector<CellSet4d>& all_current_cells,
		std::vector<CellSet4d>& refreshed_cells) {
		VEC_X_RANGE_VALIDATION(calc_index_max);

		__m256d main_cond_mask = { 0 };
		//__m256d me_val = { 0 };

		for (int i = calc_index_min / 4; i < calc_index_max / 4; i++) {
			//me_val = { 0 };
			main_cond_mask = all_current_cells[i].state.getORMask(ALIVE, FIX, MUSUME);
			if (_mm256_testz_pd(main_cond_mask, main_cond_mask) == 0) {
				for (int j = 0; j < C::max_cell_num; j++) {
					if (all_current_cells[j].hasState(ALIVE)) {
						refreshed_cells[i].w[j]= _mm256_add_pd(all_current_cells[i].w[j],
							_mm256_and_pd(
								main_cond_mask,
									_mm256_and_pd(all_current_cells[j].state.smask[ALIVE],
										Fw(all_current_cells[i].w[j], all_current_cells[i].c, all_current_cells[j].c)
									)
								)
							);
					}
				}
				
			}
			

		}
	}
}
