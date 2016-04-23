#include <math.h>
#include "define.h"
//?t?@?C??????????????????
int EPI::C::NMEMB = -1;
int EPI::C::NDER = -1;
int EPI::C::CELL_START_IDX = NMEMB + NDER; //tmp def
int EPI::Field_Data::current_cell_num = 0;
std::vector<EPI::CellSet4d> EPI::Field_Data::cells(C::max_cell_num/4);
std::vector<EPI::CellSet4d> EPI::Field_Data::next_cells(C::max_cell_num/4);
//int EPI::current_cell_num = 0;

Vec4db EPI::C::v4d_false =Vec4db(false);
Vec4db EPI::C::v4d_true = Vec4db(true);
VEC_INIT_no_inv(EPI::C, zero);
VEC_INIT(EPI::C, half);
VEC_INIT(EPI::C, one);
VEC_INIT_no_inv(EPI::C, neg_zero);
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
SET4i_d(EPI::C,NX);SET4i_d(EPI::C,NY);SET4i_d(EPI::C,NZ);
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
VEC_INIT_no_inv(EPI::Ca2P, iage_kitei);
VEC_INIT(EPI::Ca2P, Cout);


Vec4d tanh_alt(const Vec4d&x) {
	//inefficient?
    return 0.5*(1.0+tanh(x));
    //return _mm256_mul_pd(EPI::C::half_4d, _mm256_add_pd(EPI::C::one_4d, tanh_avx(x)));
}
Vec4d m256dintmod(const Vec4d &_num, const Vec4d& _den) {
    __m256d num=_num;__m256d den=_den;
    return Vec4d(_mm256_round_pd(
		_mm256_sub_pd(num,
			_mm256_mul_pd(
				_mm256_floor_pd(
					_mm256_div_pd(num, den)
					), den
				)
			),
        _MM_FROUND_NINT));
	//Round(num-Floor(num/den)*den)
}

void VSet4d::operator+=(const VSet4d &a) {
    x+=a.x;
    y+=a.y;
    z+=a.z;
}

Vec4d calc_avg8(const VSet4i& lat, const _3DScalar4d& _3DVal_4d) {
	//平均をとる番号
    std::vector<VSet4i> avg_lat(8);
    alignas(32) double avg_val[4];
    Vec4d avg;
    Vec4i avg_vec_index_x=0,avg_elem_index_x=0;
	avg_lat[0].x = lat.x;
	avg_lat[0].y = lat.y;
	avg_lat[0].z = lat.z;

	//boundary is not considered
    avg_lat[1].x = lat.x+1;
	avg_lat[1].y = lat.y;
	avg_lat[1].z = lat.z;

	avg_lat[2].x = lat.x;
    avg_lat[2].y = lat.y+1;
	avg_lat[2].z = lat.z;

	avg_lat[3].x = lat.x;
	avg_lat[3].y = lat.y;
    avg_lat[3].z = lat.z+1;

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
        avg_vec_index_x=avg_lat[k].x/const_int(4); //index of 4-pair data
        avg_elem_index_x=avg_lat[k].x-avg_vec_index_x*4; //index of the elem in 4-pair
			   //仕方ないので４つバラバラにして取得
		for (int cell_count = 0; cell_count < 4; cell_count++) {
            avg_val[cell_count]=(
                        _3DVal_4d[avg_vec_index_x[cell_count]][avg_lat[k].y[cell_count]][avg_lat[k].z[cell_count]] //vectorset
                    )[avg_elem_index_x[cell_count]];
		}

        avg +=Vec4d().load(avg_val);
	}
    return avg*0.125; // div by 8
}


namespace EPI {

	template<typename T, typename... U>
    Vec4db PState4d::getORMask(T first, U... rest) const {
        return smask[first] || getORMask(rest...);
	}
    Vec4db PState4d::getORMask() const {
        return C::v4d_false;
	}

	template<typename T, typename... U>
    Vec4db PState4d::getANDMask(T first, U... rest) const {
        return smask[first] && getANDMask(rest...);
	}
    Vec4db PState4d::getANDMask() const {
        return C::v4d_true;
	}


    void CellSet4d::get_lattice(VSet4i& out) const {
        Vec4i raw_x = round_to_int(pos.x/C::dx_4d); //inv kakero
        Vec4i raw_y = round_to_int(pos.y/C::dy_4d);
        Vec4i raw_z = round_to_int(pos.z/C::dz_4d);
        out.x = raw_x-(raw_x/const_int(C::NX))*C::NX;
        out.y = raw_y-(raw_y/const_int(C::NY))*C::NY;
        out.z = raw_z-(raw_z/const_int(C::NZ))*C::NZ;
	}
	template<typename T, typename... U>
	bool CellSet4d::hasState(T s, U... rest) const {
        return (horizontal_or(this->state.smask[s]) || hasState(rest...));
	}

	bool CellSet4d::hasState()const {
		return false;
	}
   Vec4d Ca2P::G(const VSet4d& x_vec, const VSet4d& xi_vec, const Vec4d& dcdt) {
        return select(square(x_vec.x-xi_vec.x)+square(x_vec.y-xi_vec.y)+square(x_vec.z-xi_vec.z)<=_rSq_4d && dcdt>=0,Kac_4d,0);
	}

    Vec4d Ca2P::Fc(const CellSet4d& cset, const Vec4d& B) {
		//??????????C??????????????
        Vec4d vBSq = square(B); //__m256d BSq = _mm256_mul_pd(B, B);

        Vec4d vden_Km_p = Kmu_4d+cset.P;
        Vec4d vden_K1_c = K1_4d+cset.c;
        Vec4d vden_Kg_c = Kg_4d+cset.c;//__m256d den_Kg_c = _mm256_add_pd(Kg_4d, cset.c);
        Vec4d vden_Hb_BSq = Hb_4d+vBSq;//__m256d den_Hb_BSq = _mm256_add_pd(Hb_4d, BSq);

        Vec4d vnum_mu1p = mu1_4d*cset.P;//__m256d num_mu1p = _mm256_mul_pd(mu1_4d, cset.P);
        Vec4d vnum_salphac = sub1Alpha0_4d*cset.c;//__m256d num_salphac = _mm256_mul_pd(sub1Alpha0_4d, cset.c);
        Vec4d vnum_KbcBSq = Kbc_4d*Cout_4d*vBSq;//__m256d num_KbcBSq = _mm256_mul_pd(_mm256_mul_pd(Kbc_4d,Cout_4d), BSq); //why cout?
        Vec4d vnum_gammac = gamma_4d*cset.c;//__m256d num_gammac = _mm256_mul_pd(gamma_4d, cset.c);

		//fmadd available if AVX2
        Vec4d vnum_mu0Kmu_p = vnum_mu1p+mu0_4d*vden_Km_p;//__m256d num_mu0Kmu_p = _mm256_add_pd(_mm256_mul_pd(mu0_4d, den_Km_p), num_mu1p);
        Vec4d vnum_a0K1_c_salc = vnum_salphac+alpha0_4d*vden_K1_c;//__m256d num_a0K1_c_salc = _mm256_add_pd(_mm256_mul_pd(alpha0_4d, den_K1_c), num_salphac);
        Vec4d vnum_KbcBSqKg_c__gcHB_BSq = vnum_KbcBSq*vden_Kg_c-vnum_gammac*vden_Hb_BSq;

                //__m256d num_KbcBSqKg_c__gcHb_BSq = _mm256_sub_pd(_mm256_mul_pd(num_KbcBSq, den_Kg_c), _mm256_mul_pd(num_gammac, den_Hb_BSq));

        Vec4d vden_p1 = vden_Km_p*vden_K1_c;//__m256d den_p1 = _mm256_mul_pd(den_Km_p, den_K1_c);
        Vec4d vden_p2 = vden_Kg_c*vden_Hb_BSq;//__m256d den_p2 = _mm256_mul_pd(den_Kg_c, den_Hb_BSq);
        Vec4d vnum_p1 = vden_p2*Kf_4d*cset.h*vnum_mu0Kmu_p*vnum_a0K1_c_salc;

        //__m256d num_p1 = _mm256_mul_pd(_mm256_mul_pd(den_p2, _mm256_mul_pd(Kf_4d, cset.h)), _mm256_mul_pd(num_mu0Kmu_p, num_a0K1_c_salc));
        Vec4d vnum_p2 = vden_p1*vnum_KbcBSqKg_c__gcHB_BSq;
        //__m256d num_p2 = _mm256_mul_pd(den_p1, num_KbcBSqKg_c__gcHb_BSq);
        return beta_4d+((vnum_p1*vnum_p2)/(vden_p1*vden_p2));
        //add(sub):8 mul: 17 div:1
	}
	//inline?
    Vec4d Ca2P::FP(const CellSet4d& cset, const Vec4d& A) {
        return (Kpa(cset)*A)/(H0_4d+A);
	}
    Vec4d Ca2P::Fh(const CellSet4d& cset) {
        Vec4d K2Sq = square(K2_4d);
        return (K2Sq/(cset.c+K2Sq))-cset.h;
	}

    Vec4d Ca2P::Fw(const Vec4d &wij, const Vec4d &ci, const Vec4d &cj) {
        return tanh_alt(
                    (wd0_4d-abs(ci-cj))*INV_eps_w0_4d
                    )-wij;
	}
    Vec4d Ca2P::tau_h(const CellSet4d& cset) {

        Vec4d age = select(cset.state.smask[FIX] || cset.state.smask[MUSUME],cset.ageb,cset.agek);
        Vec4d alive_val = tau_g_4d+(tau_s_4d-tau_g_4d)*tanh_alt((Ss_4d-age)*INV_delta_tau_4d);
        return select(cset.state.smask[ALIVE],alive_val,select(cset.state.smask[FIX] || cset.state.smask[MUSUME],tau_s_4d,zero_4d));
	}
    Vec4d Ca2P::In(const CellSet4d& cset) {
        return select(
                    cset.state.smask[FIX] || cset.state.smask[MUSUME],
                    iage_kitei_4d,
                    tanh_alt((cset.agek-Ss_4d)*INV_delta_I_4d)
                    );
	}

    Vec4d Ca2P::Kpa(const CellSet4d& cset) {
        Vec4d age = select(cset.state.smask[FIX] || cset.state.smask[MUSUME],cset.ageb,cset.agek);
        Vec4d alive_val = kg_4d+(ks_4d-kg_4d)*tanh_alt((Ss_4d-age)*INV_delta_k_4d);
        return select(cset.state.smask[ALIVE],alive_val,select(cset.state.smask[FIX] || cset.state.smask[MUSUME],ks_4d,zero_4d)); //zero->exception on somewhere
	}

	void Ca2P::refresh_Ca(const CUBE& calc_area,
		const _3DScalar4d& currentCa,
		const std::vector<CellSet4d>& all_cells,
		_3DScalar4d& nextCa) {
		//x is compressed by 4
		VEC_X_RANGE_VALIDATION(calc_area.x.min); VEC_X_RANGE_VALIDATION(calc_area.x.max);

		int prev_x_idx, next_x_idx, prev_y_idx, next_y_idx, prev_z_idx, next_z_idx;

		double _i, _i1, _i2, _i3, _i4;
        Vec4d _sum = 0;
        Vec4d x_next_sub, x_prev_sub, y_next_sub, y_prev_sub, z_next_sub, z_prev_sub, sub_sum;
        Vec4d x_tmp;
		VSet4d CaPos;
		//currently boundary conditon is not considered
		for (int i = calc_area.x.min / 4; i <= calc_area.x.max / 4; i++) {
			_i = i * 4;
			_i1 = _i*dx; _i2 = _i1 + dx; _i3 = _i2 + dx; _i4 = _i3 + dx;
            CaPos.x = Vec4d(_i1, _i2, _i3, _i4);
			//periodic boundary condition for x and y
			prev_x_idx = ((i - 4 + NX) % NX) / 4;
			next_x_idx = ((i + 4 + NX) % NX) / 4;
			for (int j = calc_area.y.min; j <= calc_area.y.max; j++) {
                CaPos.y = Vec4d(i*dy);
				prev_y_idx = (j - 1 + NY) % NY;
				next_y_idx = (j + 1 + NY) % NY;
				for (int k = calc_area.z.min; k <= calc_area.z.max; k++) {
                    _sum=0;
                    CaPos.z = Vec4d(i*dz);
					//neumann
					prev_z_idx = k == 0 ? 1 : k - 1;
					next_z_idx = k == NZ ? NZ - 1 : k + 1;
					for (int l = CELL_START_IDX; l < Field_Data::current_cell_num; l++) {
                        _sum += G(CaPos, all_cells[l].pos, all_cells[l].diff_c);
					}
					//there is a more efficient way
                    x_next_sub = Vec4d(currentCa[i][j][k] [1], currentCa[i][j][k] [2], currentCa[i][j][k] [3], currentCa[next_x_idx][j][k] [0])-currentCa[i][j][k];
                    x_prev_sub = Vec4d(currentCa[prev_x_idx][j][k] [3], currentCa[i][j][k] [0], currentCa[i][j][k] [1], currentCa[i][j][k] [2])-currentCa[i][j][k];
                    y_next_sub = currentCa[i][next_y_idx][k]-currentCa[i][j][k];
                    y_prev_sub = currentCa[i][prev_y_idx][k] -currentCa[i][j][k];
                    z_next_sub = currentCa[i][j][next_z_idx]-currentCa[i][j][k];
                    z_prev_sub = currentCa[i][j][prev_z_idx]- currentCa[i][j][k];
                    sub_sum=x_next_sub+x_prev_sub+y_next_sub+y_prev_sub+z_next_sub+z_prev_sub;

                    nextCa[i][j][k] = currentCa[i][j][k]+dt_ca_4d*((dA_4d*sub_sum/dxSq_4d)+_sum-Kaa_4d*currentCa[i][j][k]);


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

        VSet4i lat;
        Vec4d neg_me_dead_val = 0, react_alive_val = 0, me_val = 0 , react_val = 0 ;//dt unmultiplied
        Vec4d first_calc_val =  0 , second_calc_val =  0 ;
        Vec4db react_mask =  false ;
        Vec4d In_val =  0 ;
        Vec4db main_cond_mask =  false ;
		for (int i = calc_index_min / 4; i < calc_index_max / 4; i++) {
            neg_me_dead_val =  0 ;
            react_alive_val =  0 ;
            me_val =  0 ;
            react_val =  0 ;

			//first calc
			if (all_current_cells[i].hasState(DEAD)) {
				//Kpp*Pi
                neg_me_dead_val = Kpp_4d*all_current_cells[i].P;
				for (int j = 0; j < C::max_cell_num / 4; j++) {
					if (all_current_cells[j].hasState(ALIVE)) { //これでは全ての細胞に対して計算する
						//やっぱり重いなら変える
						//Pj-Pi
                        react_alive_val += select(all_current_cells[j].state.smask[ALIVE],all_current_cells[j].P- all_current_cells[i].P,0);
					}
				}
				//dp*Sum(Pj-Pi)
                react_alive_val *= dP_4d;
			}
            first_calc_val = select(all_current_cells[i].state.smask[DEAD],react_alive_val- neg_me_dead_val,0);
			//先にDEADマスクをかけずにneg_me_dead_val+react_alive_valを計算してからかけたほうが良い?
			//->変わらない,けど見やすいかも

            //better:use ALIVE||FIX||MUSUME MASK for check cond and later masking
            main_cond_mask = all_current_cells[i].state.getORMask(ALIVE,FIX,MUSUME);

            if(horizontal_or(main_cond_mask)){
                all_current_cells[i].get_lattice(lat);
                me_val = FP(all_current_cells[i],calc_avg8(lat,currentCa))-Kpp_4d*all_current_cells[i].P;
                In_val=In(all_current_cells[i]);
                for(int j=0;j<C::max_cell_num/4;j++){
                    react_mask=all_current_cells[j].state.getORMask(ALIVE,DEAD,FIX,MUSUME);
                    if(horizontal_or(react_mask)){
                        react_val +=select(react_mask,In_val*all_current_cells[i].w[j]*(all_current_cells[j].P-all_current_cells[i].P),0);
                    }
                }
                react_val*=dP_4d;
            }

            second_calc_val = select(main_cond_mask,me_val+react_val,0);
            refreshed_cells[i].P=all_current_cells[i].P+C::dt_ca_4d*(first_calc_val+second_calc_val);
		}
	}

	void Ca2P::refresh_c_i(int calc_index_min, int calc_index_max,
		const _3DScalar4d& currentB,
		const std::vector<CellSet4d>& all_current_cells,
        std::vector<Vec4d>& diff_c_out,
		std::vector<CellSet4d>& refreshed_cells) {
		VEC_X_RANGE_VALIDATION(calc_index_min);
		VEC_X_RANGE_VALIDATION(calc_index_max);
        Vec4db main_cond_mask =  false ;
        Vec4db react_mask =  false ;
        Vec4d In_val =  0 , react_val =  0 ;
        VSet4i lat;
        std::fill(diff_c_out.begin(), diff_c_out.end(), 0);
		for (int i = calc_index_min / 4; i < calc_index_max / 4; i++) {
            react_val =  0 ;
            main_cond_mask =all_current_cells[i].state.getORMask(ALIVE,FIX,MUSUME);

            if (horizontal_or(main_cond_mask)) {
				all_current_cells[i].get_lattice(lat);
				//default diffu=0
				In_val = In(all_current_cells[i]);
				for (int j = 0; j < C::max_cell_num / 4; j++) {
                    react_mask=all_current_cells[j].state.getORMask(ALIVE,DEAD,FIX,MUSUME);
                    if (horizontal_or(react_mask)) {
                        react_val+=select(react_mask,In_val*all_current_cells[i].w[j]*(all_current_cells[j].c-all_current_cells[i].c),0);
					}
                }
				//mask over all value,
                diff_c_out[i] = select(main_cond_mask,Fc(all_current_cells[i], calc_avg8(lat, currentB))+dc_4d+react_val,0);
			}
			//0 when the cond not satisfied
			//see -> u_ave[j]=u[j]=0.0;
            refreshed_cells[i].c = select(main_cond_mask,all_current_cells[i].c+dt_ca_4d*diff_c_out[i],0);
		}
	}

	void Ca2P::refresh_h_i(int calc_index_min, int calc_index_max,
		const std::vector<CellSet4d>& all_current_cells,
		std::vector<CellSet4d>& refreshed_cells) {
		VEC_X_RANGE_VALIDATION(calc_index_min);
		VEC_X_RANGE_VALIDATION(calc_index_max);

        Vec4db main_cond_mask =  false ;
        Vec4d me_val =  0 ;
		for (int i = calc_index_min / 4; i < calc_index_max / 4; i++) {
            me_val = 0;
			main_cond_mask = all_current_cells[i].state.getORMask(ALIVE, FIX, MUSUME);
            if (horizontal_or(main_cond_mask)) {
                me_val += select(main_cond_mask,Fh(all_current_cells[i])/tau_h(all_current_cells[i]),0);
			}
            refreshed_cells[i].h = all_current_cells[i].h+select(main_cond_mask,dt_ca_4d*me_val,0);

		}
	}

	void Ca2P::refresh_w_i_j(int calc_index_min, int calc_index_max,
		const std::vector<CellSet4d>& all_current_cells,
		std::vector<CellSet4d>& refreshed_cells) {
		VEC_X_RANGE_VALIDATION(calc_index_max);

        Vec4db main_cond_mask = false;
		//__m256d me_val = { 0 };

		for (int i = calc_index_min / 4; i < calc_index_max / 4; i++) {
			//me_val = { 0 };
			main_cond_mask = all_current_cells[i].state.getORMask(ALIVE, FIX, MUSUME);
            if (horizontal_or(main_cond_mask)) {
				for (int j = 0; j < C::max_cell_num; j++) {
					if (all_current_cells[j].hasState(ALIVE)) {

                        refreshed_cells[i].w[j]=all_current_cells[i].w[j]+
                                select(main_cond_mask && all_current_cells[j].state.smask[ALIVE],
                                       Fw(all_current_cells[i].w[j], all_current_cells[i].c, all_current_cells[j].c),
                                       0);
					}
				}
				
			}
			

		}
	}
}
