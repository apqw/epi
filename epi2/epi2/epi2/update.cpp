//#include <stdlib.h>
#include <stdio.h>
#include "update.h"

#include <omp.h>
#include <cassert>
#include <vector>
#include "SFMT.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <cstring>

#define max(a,b) ((a)<=(b)?(b):(a))
#define min(a,b) ((a)<=(b)?(a):(b))

#define ret_x(x1,x2) (ret__((x1),(x2),CONST::Lx))
#define ret_y(y1,y2) (ret__((y1),(y2),CONST::Ly))

#define TEST

//http://www.geometrictools.com/Documentation/FastInverseSqrt.pdf
//dont use this in production 
inline double t_sqrtD(const double &x)
{
	double         xHalf = 0.5 * x;
	long long int  tmp = 0x5FE6EB50C7B537AAl - (*(long long int*)&x >> 1);
	double         xRes = *(double*)&tmp;

	xRes *= (1.5 - (xHalf * xRes * xRes));
	xRes *= ( 1.5 - ( xHalf * xRes * xRes ) );//more precise
	return xRes * x;
}


inline double ret__(double x1, double x2, double axis_len)
/* choose one of LX-periodic values of x2 such that
|x1-x2| takes minimum
*/
{
	double min, tmp;
	min = fabs(x2 - x1);
	tmp = x2;
	if (min > fabs(x2 - axis_len - x1)) {
		min = fabs(x2 - axis_len - x1);
		tmp = x2 - axis_len;
	}
	if (min > fabs(x2 + axis_len - x1)) {
		min = fabs(x2 + axis_len - x1);
		tmp = x2 + axis_len;
	}
	return tmp;
}



inline double calc_avg8(const double arr[CONST::NX+1][CONST::NY+1][CONST::NZ+1], int ix, int iy, int iz) {
	return 0.125*(arr[ix][iy][iz] + arr[ix][iy][iz + 1] +
		+arr[ix][iy + 1][iz] + arr[ix][iy + 1][iz + 1]
		+ arr[ix + 1][iy][iz] + arr[ix + 1][iy][iz + 1]
		+ arr[ix + 1][iy + 1][iz] + arr[ix + 1][iy + 1][iz + 1]);
}

inline double fu(double u, double v, double p, double B) {
	//kgとKgがある　何故
	return CONST::Ca2P::Kf*(CONST::Ca2P::mu0 + CONST::Ca2P::mu1*p / (p + CONST::Ca2P::Kmu))*(CONST::Ca2P::alpha0 + CONST::Ca2P::sub1_alpha0*u / (CONST::Ca2P::K1 + u))*v
		- CONST::Ca2P::gamma*u / (CONST::Ca2P::Kg + u) + CONST::Ca2P::beta + CONST::Ca2P::Kbc*B*B*CONST::Ca2P::Cout / (CONST::Ca2P::Hb + B*B);
}

inline double fh(double u, double v) {
	return (CONST::Ca2P::K2Sq / (CONST::Ca2P::K2Sq + u*u)) - v;
}

inline double tau_h(double age, CELL_STATE state) {
	if (state == ALIVE) {
		return CONST::Ca2P::tau_g + 0.5*(CONST::Ca2P::tau_s - CONST::Ca2P::tau_g)*(1 + tanh((CONST::THRESH_SP - age) / CONST::Ca2P::delta_tau));
	}
	else if (state == FIX || state == MUSUME) {
		return CONST::Ca2P::tau_s;
	}
	else {
		throw "invalid value in fv";
		exit(1);
	}
}

inline double Kpa(double age, CELL_STATE state) {
	if (state == ALIVE) {
		return CONST::Ca2P::kg + 0.5*(CONST::Ca2P::ks - CONST::Ca2P::kg)*(1 + tanh((CONST::THRESH_SP - age) / CONST::Ca2P::delta_k));
	}
	else if (state == FIX || state == MUSUME) {
		return CONST::Ca2P::ks;
	}
	else {
		return 0;
	}
}

inline double IAG(double age, CELL_STATE state) {
	if (state == FIX || state == MUSUME) {
		return CONST::Ca2P::iage_kitei;
	}
	else {
		return 0.5*(1.0 + tanh((age - CONST::THRESH_SP) / CONST::Ca2P::delta_I));
	}
}

inline double fp(double a, double age, CELL_STATE state) {
	return Kpa(age, state)*a / (CONST::Ca2P::H0 + a);
}

inline double fw(double diff, double w) {
	return 0.5*(1.0 + tanh((CONST::Ca2P::wd0 - diff) / CONST::Ca2P::eps_w0)) - w;
}
void update_cell_internal(const ENV& env, ENV& updated_env) {
	//

    //for (int i = 0; i < CONST::Ca2P::ca_N; i++) {

        double *diffu = updated_env.cell_diff_c;
        //printf("%d\n", i);
#pragma omp parallel for
		for (int j = env.NMEMB+env.NDER; j < env.current_cell_num; j++) {

			//calc p (pre)
			double p_diff = 0;
			if (env.cell_state[j] == DEAD) {
				p_diff -= CONST::Ca2P::Kpp*env.cell_P[j];
				double unmult_dp = 0;
				int connect_alive_count = 0;
				for (int k = 0; k < env.cell_connected_num[j]; k++) {
					int connected_cell_index = env.cell_connected_index[j][k];
					if (connected_cell_index > -1 && env.cell_state[connected_cell_index] == ALIVE) {
						connect_alive_count++;
						unmult_dp += env.cell_P[connected_cell_index];
					}
				}
				p_diff += CONST::Ca2P::dP*(unmult_dp - env.cell_P[j] * connect_alive_count);
			}
			
			double v_diff = 0;
			if (env.cell_state[j] & CONST::ALIVE_FIX_MUSUME) {
				int ix = env.cell_x_index[j], iy = env.cell_y_index[j], iz = env.cell_z_index[j];
				double avg_Ca = calc_avg8(env.Ca2P_value, ix, iy, iz);
				double avg_B = calc_avg8(env.B_value, ix, iy, iz);
				double age = (env.cell_state[j] &(FIX|MUSUME)) ? env.cell_ageb[j] : env.cell_agek[j];
				diffu[j] = fu(env.cell_c[j], env.cell_h[j], env.cell_P[j], avg_B);
				v_diff += fh(env.cell_c[j], env.cell_h[j]) / tau_h(age, env.cell_state[j]);
				p_diff += fp(avg_Ca, age, env.cell_state[j])-CONST::Ca2P::Kpp*env.cell_P[j];
			}
			double diffu_unmult = 0;
			double IAG_val = IAG(env.cell_agek[j], env.cell_state[j]);
			double p_diff_unmult = 0;
            //double w_diff = 0;
			for (int k = 0; k < env.cell_connected_num[j]; k++) {
				int connected_cell_index = env.cell_connected_index[j][k];
				if (env.cell_state[connected_cell_index] & (CONST::ALIVE_FIX_MUSUME | DEAD)) {
					double u_tmp_diff = env.cell_c[connected_cell_index] - env.cell_c[j];
					diffu_unmult += env.cell_w[j][k] * u_tmp_diff;
					p_diff_unmult+= env.cell_w[j][k] * (env.cell_P[connected_cell_index] - env.cell_P[j]);

					if (env.cell_state[connected_cell_index] == ALIVE) {
						updated_env.cell_w[j][k] = env.cell_w[j][k]+ CONST::dt_ca*fw(fabs(u_tmp_diff), env.cell_w[j][k]);
					}
				}
			}
			diffu[j] += CONST::Ca2P::dc*IAG_val*diffu_unmult;
			updated_env.cell_c[j] = env.cell_c[j] + CONST::dt_ca*diffu[j];
            updated_env.cell_ca_ave[j] +=updated_env.cell_c[j];
			updated_env.cell_P[j] = env.cell_P[j] + CONST::dt_ca*(p_diff + CONST::Ca2P::dP*IAG_val*p_diff_unmult);
			updated_env.cell_h[j] = env.cell_h[j] + CONST::dt_ca*v_diff;
			

		}


    //}
}
#define in_loop_z do{\
map = current_env->cell_map[j][k][l];\
mapexist = isNOTMINUS(map);\
map = mapexist*map;\
flg = current_env->cell_state[map] & CONST::ALIVE_FIX_MUSUME;\
result = isNOTZERO(flg)&mapexist;\
diffu = result*current_env->cell_diff_c[map*result];\
\
centerCa2P = current_env->Ca2P_value[j][k][l];\
\
astim = current_env->air_stim[j][k][l];\
mul1 = _mm256_set_pd(\
	current_env->cell_map2[prev_x][k][l],\
	current_env->cell_map2[next_x][k][l],\
	current_env->cell_map2[j][prev_y][l],\
	current_env->cell_map2[j][next_y][l]);\
mul2 = _mm256_set_pd(\
	current_env->cell_map2[j][k][prev_z],\
	current_env->cell_map2[j][k][next_z],\
	0,\
	0\
	);\
sub1_1 = _mm256_set_pd(\
	current_env->Ca2P_value[prev_x][k][l],\
	current_env->Ca2P_value[next_x][k][l],\
	current_env->Ca2P_value[j][prev_y][l],\
	current_env->Ca2P_value[j][next_y][l]\
	);\
sub1_2 = _mm256_set_pd(\
	current_env->Ca2P_value[j][k][prev_z],\
	current_env->Ca2P_value[j][k][next_z],\
	0,\
	0\
	);\
sub2 = _mm256_set1_pd(centerCa2P);\
\
_mm256_store_pd(diffuse_coef, _mm256_add_pd(_mm256_mul_pd(mul1, _mm256_sub_pd(sub1_1, sub2)), _mm256_mul_pd(mul2, _mm256_sub_pd(sub1_2, sub2))));\
updated_env->Ca2P_value[j][k][l] = centerCa2P\
+ CONST::dt_ca*(CONST::Ca2P::dA*(\
	diffuse_coef[0] + diffuse_coef[1] + diffuse_coef[2] + diffuse_coef[3])\
	*CONST::inv_dxSq + CONST::Ca2P::Kac*diffu - CONST::Ca2P::Kaa*centerCa2P + astim);\
}while(0)
			
//storeが遅い

void update_ca(const ENV *current_env,ENV *updated_env){
    int iz_bound=(int)((current_env->zzmax+CONST::FAC_MAP*CONST::R_max)*CONST::inv_dz);
#ifdef TEST
	//iz_bound = 20;
#endif
	//register double out[CONST::NX][CONST::NY][CONST::NZ];
	//printf("iz_bound:%d\n", iz_bound);
	auto out = updated_env->Ca2P_value;
#pragma omp parallel for
		for (int j = 0; j < CONST::NX; j++) {
			__m256d mul1, mul2, sub1_1, sub1_2, sub2;
			alignas(32) double diffuse_coef[4];
			double xx;
			double diffu, centerCa2P, astim;
			int prev_z, next_z, prev_x, next_x, prev_y, next_y, map;
			unsigned int mapexist, flg, result;

			prev_x = (j - 1 + CONST::NX) % CONST::NX;
			next_x = (j + 1 + CONST::NX) % CONST::NX;
			//double xxa[CONST::NX][CONST::NY][CONST::NZ];
			for (int k = 0; k < CONST::NY; k++) {
				prev_y = (k - 1 + CONST::NY) % CONST::NY;
				next_y = (k + 1 + CONST::NY) % CONST::NY;

				//double *zza = yya[k];
				for (int l = 1; l < iz_bound; l++) {
					//prev_z = l==0?1:l - 1;
					//next_z = l==CONST::NZ?CONST::NZ-1:l+1;
					//next_z = l + 1;
					//test
					prev_z = l - 1; next_z = l + 1;
					//in_loop_z;
					map = current_env->cell_map[j][k][l];
					mapexist = isNOTMINUS(map);
					map = mapexist*map;
					flg = current_env->cell_state[map] & CONST::ALIVE_FIX_MUSUME;
					result = isNOTZERO(flg)&mapexist;
					diffu = result*current_env->cell_diff_c[map*result];

					centerCa2P = current_env->Ca2P_value[j][k][l];

					astim = current_env->air_stim[j][k][l];
					//slow I/O
					out[j][k][l] = centerCa2P
						+ CONST::dt_ca*(CONST::Ca2P::dA*(
							current_env->cell_map2[prev_x][k][l] * (current_env->Ca2P_value[prev_x][k][l] - centerCa2P)
							+ current_env->cell_map2[next_x][k][l] * (current_env->Ca2P_value[next_x][k][l] - centerCa2P)
							+ current_env->cell_map2[j][prev_y][l] * (current_env->Ca2P_value[j][prev_y][l] - centerCa2P)
							+ current_env->cell_map2[j][next_y][l] * (current_env->Ca2P_value[j][next_y][l] - centerCa2P)
							+ current_env->cell_map2[j][k][prev_z] * (current_env->Ca2P_value[j][k][prev_z] - centerCa2P)
							+ current_env->cell_map2[j][k][next_z] * (current_env->Ca2P_value[j][k][next_z] - centerCa2P)
							)*CONST::inv_dxSq + CONST::Ca2P::Kac*diffu - CONST::Ca2P::Kaa*centerCa2P+astim);


				}
				//memcpy(yya[k], zza, iz_bound*sizeof(double));
			}
			//memcpy
		}
    //periodic
    for(int l=0;l<=iz_bound;l++){
        for(int j=0;j<CONST::NX;j++)updated_env->Ca2P_value[j][CONST::NY][l] = updated_env->Ca2P_value[j][0][l];
        for(int j=0;j<=CONST::NY;j++)updated_env->Ca2P_value[CONST::NX][j][l] = updated_env->Ca2P_value[0][j][l];
    }
}

void update_map(ENV *update_env){
    //optimized
	std::fill(&(update_env->cell_map[0][0][0]),&(update_env->cell_map[0][0][0])+CONST::NX*CONST::NY*CONST::NZ,-1);
    std::fill(&(update_env->cell_map2[0][0][0]),&(update_env->cell_map2[0][0][0])+CONST::NX*CONST::NY*CONST::NZ,0);
//大した量ではない？    
#pragma omp parallel for
    for(int i=0;i<update_env->current_cell_num;i++){
        if(update_env->cell_state[i]!=DISA){
            int ix = update_env->cell_x_index[i];
            int iy = update_env->cell_y_index[i];
            int iz = update_env->cell_z_index[i];
            double rad=update_env->cell_radius[i];
            double radSq=rad*rad;
            int irx = (int)(rad*CONST::inv_dx);
            int iry = (int)(rad*CONST::inv_dy);
            int irz = (int)(rad*CONST::inv_dz);
            int izmin = max(1,iz-irz)-iz;int izmax=min(CONST::NZ-1,iz+irz)-iz;
            int iymin=-iry;int iymax=iry;
            int ixmin=irx;int ixmax=irx;

            for(int j=izmin;j<=izmax;j++){ //width! not coord
                double zradSq=j*j*CONST::dzSq;
				int ipz= iz + j;
                for(int k=iymin;k<=iymax;k++){
                    double yradSq=k*k*CONST::dySq;
                    int imy=k+iy;
                    int ipy=imy<0?imy+CONST::NY:
                                  imy>=CONST::NY?
                                               imy-CONST::NY:
                                                imy;//when imx<-NX or imx>2*NX ,this part will be broken
                    for(int l=ixmin;l<=ixmax;l++){
                        int imx=l+ix;
                        int ipx=imx<0?imx+CONST::NX:
                                      imx>=CONST::NX?
                                                   imx-CONST::NX:
                                                    imx;
                        if(l*l*CONST::dxSq+yradSq+zradSq<=radSq){
                            update_env->cell_map[ipx][ipy][ipz]=i;//slow cache?
                        }
                    }
                }
            }
            if(update_env->cell_state[i] & CONST::ALIVE_FIX_MUSUME){
                rad=2.0*update_env->cell_radius[i];
                radSq=rad*rad;
              irx = (int)(rad*CONST::inv_dx);
                iry = (int)(rad*CONST::inv_dy);
                irz = (int)(rad*CONST::inv_dz);
                izmin = max(1,iz-irz)-iz;izmax=min(CONST::NZ-1,iz+irz)-iz;
                iymin=-iry;iymax=iry;
                ixmin=irx;ixmax=irx;


                for(int j=izmin;j<=izmax;j++){ //width! not coord
                    double zradSq=j*j*CONST::dzSq;
					int ipz = iz + j;
                    for(int k=iymin;k<=iymax;k++){
                        double yradSq=k*k*CONST::dySq;
                        int imy=k+iy;
                        int ipy=imy<0?imy+CONST::NY:
                                      imy>=CONST::NY?
                                                   imy-CONST::NY:
                                                    imy;//when imx<-NX or imx>2*NX ,this part will be broken
                        for(int l=ixmin;l<=ixmax;l++){
                            int imx=l+ix;
                            int ipx=imx<0?imx+CONST::NX:
                                          imx>=CONST::NX?
                                                       imx-CONST::NX:
                                                        imx;
                            if(l*l*CONST::dxSq+yradSq+zradSq<=radSq){
                                update_env->cell_map2[ipx][ipy][ipz]=1; //slow cache?
                            }
                        }
                    }
                }
            }
        }
    }

    for(int l=0;l<CONST::NZ;l++){
        for(int j=0;j<CONST::NX;j++){
            update_env->cell_map[j][CONST::NY][l]=update_env->cell_map[j][0][l];
            update_env->cell_map2[j][CONST::NY][l]=update_env->cell_map2[j][0][l];
        }
        for(int k=0;k<=CONST::NY;k++){
            update_env->cell_map[CONST::NX][k][l]=update_env->cell_map[0][k][l];
            update_env->cell_map2[CONST::NX][k][l]=update_env->cell_map2[0][k][l];
        }
    }
}

void update_all(ENV **current_env, ENV **updated_env){
	static int sw = 0;
    update_map(*current_env);
    //reset_ca_ave_in_cell(*updated_env); //do not use current value
   // set_air_stim(*current_env);
    ENV* old_env,*new_env,*tmp_env_ptr;
    old_env=*current_env;new_env=*updated_env;
	update_cell_dyn(old_env, new_env);
	//test
	sw++;
    if (sw > 20 && false) {
		for (int i = 0; i < CONST::Ca2P::ca_N; i++) {
			//printf("%d\n", i);
			update_cell_internal(*old_env, *new_env);
			update_ca(old_env, new_env);

			tmp_env_ptr = new_env;
			new_env = old_env;
			old_env = tmp_env_ptr; //old is latest now
		}
		for (int i = old_env->NMEMB + old_env->NDER; i < old_env->current_cell_num; i++) {
			if (old_env->cell_state[i] & CONST::ALIVE_FIX_MUSUME) {
				old_env->cell_ca_ave[i] /= CONST::Ca2P::ca_N;
			}
			else {
				old_env->cell_ca_ave[i] = old_env->cell_c[i] = 0.0;
			}
		}
		sw = 0;
	}
}

void reset_ca_ave_in_cell(ENV *update_env){
    for(int i=update_env->NMEMB+update_env->NDER;i<update_env->current_cell_num;i++){
        if(update_env->cell_state[i]&CONST::ALIVE_FIX_MUSUME){
            update_env->cell_ca_ave[i]=0.0;
        }
    }
}

void set_air_stim(ENV *current_env){
    for(int i=0;i<=CONST::NX;i++){
        for(int j=0;j<=CONST::NY;j++){
            for(int k=0;k<=CONST::NZ;k++){
                current_env->air_stim[i][j][k]=0.0;
                int jj=current_env->cell_map[i][j][k];
                if(jj!=-1){
                    for(int l=0;l<current_env->cell_connected_num[jj];l++){
                        if(current_env->cell_state[current_env->cell_connected_index[jj][l]]==AIR){
                            current_env->air_stim[i][j][k]=CONST::Ca2P::AIR_STIM;
                            break;
                        }
                    }
                }
            }
        }
    }
}

////////////////////////////////////////////////////
void set_memb_indices(ENV* env) {
	//jを与えて、曲げ弾性に使う6点を得られるような配列の作成？
	//何度も実行するわけではないので処理速度は適当で良い?
	//rr,uuも追加
	for (int j = 0; j < env->NMEMB; j++) {
		int jj = j%CONST::NMX; 
		int kk = j / CONST::NMX;
		env->memb_idx_jl[j] = j + (jj == 0 ? CONST::NMX : 0) - 1;
		env->memb_idx_jll[j] = j + (jj <= 1 ? CONST::NMX : 0) - 2;
		env->memb_idx_jr[j] = j + (jj == CONST::NMX - 1 ? CONST::NMX : 0) + 1;
		env->memb_idx_jrr[j] = j + (jj == CONST::NMX - 2 ? CONST::NMX : 0) + 2;

		env->memb_idx_jb[j] = j + (kk == 0 ? CONST::NMY : 0)*CONST::NMX - CONST::NMX;
		env->memb_idx_jbb[j] = j + (kk <= 1 ? CONST::NMY : 0)*CONST::NMX - 2 * CONST::NMX;
		env->memb_idx_ju[j] = j - (kk == CONST::NMY - 1 ? CONST::NMY : 0)*CONST::NMX + CONST::NMX;
		env->memb_idx_juu[j] = j- (kk == CONST::NMY - 2 ? CONST::NMY : 0)*CONST::NMX + 2*CONST::NMX;

		if (env->memb_idx_jl[j] < 0 || env->memb_idx_jl[j] >= env->NMEMB) { printf("error: illegal index\n"); exit(1); }
		if (env->memb_idx_jll[j] < 0 || env->memb_idx_jll[j] >= env->NMEMB) { printf("error: illegal index\n"); exit(1); }
		if (env->memb_idx_jr[j] < 0 || env->memb_idx_jr[j] >= env->NMEMB) { printf("error: illegal index\n"); exit(1); }
		if (env->memb_idx_jrr[j] < 0 || env->memb_idx_jrr[j] >= env->NMEMB) { printf("error: illegal index\n"); exit(1); }
		if (env->memb_idx_jb[j] < 0 || env->memb_idx_jb[j] >= env->NMEMB) { printf("error: illegal index\n"); exit(1); }
		if (env->memb_idx_jbb[j] < 0 || env->memb_idx_jbb[j] >= env->NMEMB) { printf("error: illegal index\n"); exit(1); }
		if (env->memb_idx_ju[j] < 0 || env->memb_idx_ju[j] >= env->NMEMB) { printf("error: illegal index\n"); exit(1); }
		if (env->memb_idx_juu[j] < 0 || env->memb_idx_juu[j] >= env->NMEMB) { printf("error: illegal index\n"); exit(1); }
	}
}

inline double ubend(const __m256d& xi, const __m256d& xj, const __m256d& xk) {
	__m256d sub1 = _mm256_sub_pd(xi, xj);
	__m256d sub2 = _mm256_sub_pd(xj, xk);
	__m256d cross = _mm256_mul_pd(sub1, sub2);
	__m256d distSq1 = _mm256_mul_pd(sub1, sub1);
	__m256d distSq2 = _mm256_mul_pd(sub2, sub2);
	alignas(32) double dSqArr1[4], dSqArr2[4],crossArr[4];
	_mm256_store_pd(crossArr,cross);
	_mm256_store_pd(dSqArr1, distSq1); _mm256_store_pd(dSqArr2, distSq2);
	double cross_val = crossArr[0] + crossArr[1] + crossArr[2];
	return 0.5*CONST::KBEND*cross_val*cross_val / ((dSqArr1[0] + dSqArr1[1] + dSqArr1[2])*(dSqArr2[0] + dSqArr2[1] + dSqArr2[2]));
}
/* return x1-x2 taking into account periodic bc */
double perio_diff(double x1, double x2, double lb)
{
	double dum;

	if (fabs(x1 - x2) < 0.5*lb) return x1 - x2;
	else {
		if (x1 > x2) {
			dum = x1 - lb;
			return dum - x2;
		}
		else {
			dum = x2 - lb;
			return x1 - dum;
		}
	}
}

inline double bend_force_sqr(double n[], double m[], double ipn[], double ipm[], double inv_dn[], double inv_dm[],
	int j, int jr, int jl, int jll, int ju, int jb, int jbb)
{
	return CONST::KBEND*(
		-(1.0 - ipn[j])*(n[jr] - ipn[j] * n[j])*inv_dn[j]
		+ (1.0 - ipn[jl])*(-(n[jl] - ipn[jl] * n[j])*inv_dn[j]
			+ (n[j] - ipn[jl] * n[jl])*inv_dn[jl])
		+ (1.0 - ipn[jll])*(n[jll] - ipn[jll] * n[jl])*inv_dn[jl]

		- (1.0 - ipm[j])*(m[ju] - ipm[j] * m[j]) *inv_dm[j]
		+ (1.0 - ipm[jb])*(-(m[jb] - ipm[jb] * m[j])*inv_dm[j]
			+ (m[j] - ipm[jb] * m[jb]) *inv_dm[jb])
		+ (1.0 - ipm[jbb])*(m[jbb] - ipm[jbb] * m[jb])*inv_dm[jb]);
}

void bend(const ENV* env,ENV* new_env,double nx[],double ny[],double nz[], double mx[], double my[], double mz[],
	double inv_dn[],double inv_dm[],double ipn[],double ipm[]) {
	double diffx, diffy, diffz,cellx,celly,cellz;
	double memb_jr_x, memb_jr_y, memb_jr_z, memb_ju_x, memb_ju_y, memb_ju_z,inv_tmp;
	for (int j = 0; j < env->NMEMB; j++) {
		cellx = env->cell_x[j];
		celly = env->cell_y[j];
		cellz = env->cell_z[j];
		memb_jr_x = env->cell_x[env->memb_idx_jr[j]];
		memb_jr_y = env->cell_y[env->memb_idx_jr[j]];
		memb_jr_z = env->cell_z[env->memb_idx_jr[j]];
		memb_ju_x = env->cell_x[env->memb_idx_ju[j]];
		memb_ju_y = env->cell_y[env->memb_idx_ju[j]];
		memb_ju_z = env->cell_z[env->memb_idx_ju[j]];
		diffx = memb_jr_x - ret_x(memb_jr_x, cellx);
		diffy = memb_jr_y - ret_y(memb_jr_y, celly);
		diffz = memb_jr_z - cellz;

		//replace this real sqrt func in production env
		inv_tmp = inv_dn[j] = t_sqrtD(diffx*diffx+diffy*diffy+diffz*diffz);

		assert(inv_tmp <= 1.0e+5);

		nx[j] = perio_diff(memb_jr_x, cellx, CONST::Lx)*inv_tmp;
		ny[j] = perio_diff(memb_jr_y, celly, CONST::Ly)*inv_tmp;
		nz[j] = diffz*inv_tmp;

		diffx = memb_ju_x - ret_x(memb_ju_x, cellx);
		diffy = memb_ju_y - ret_y(memb_ju_y, celly);
		diffz = memb_ju_z - cellz;
		inv_tmp=inv_dm[j] = t_sqrtD(diffx*diffx + diffy*diffy + diffz*diffz);

		assert(inv_tmp <= 1.0e+5);

		mx[j] = perio_diff(memb_ju_x, cellx, CONST::Lx)*inv_tmp;
		my[j] = perio_diff(memb_ju_y, cellx, CONST::Ly)*inv_tmp;
		mz[j] = diffz*inv_tmp;

	}
	for (int j = 0; j < env->NMEMB; j++) {
		ipn[j] = nx[env->memb_idx_jr[j]] * nx[j] + ny[env->memb_idx_jr[j]] * ny[j] + nz[env->memb_idx_jr[j]] * nz[j];
		ipm[j] = mx[env->memb_idx_ju[j]] * mx[j] + my[env->memb_idx_ju[j]] * my[j] + mz[env->memb_idx_ju[j]] * mz[j];
	}

	for (int j = 0; j < env->NMEMB; j++) {
		new_env->cell_x[j] =cellx+ CONST::dt_cell*bend_force_sqr(nx, mx, ipn, ipm, inv_dn, inv_dm, j, env->memb_idx_jr[j], env->memb_idx_jl[j], env->memb_idx_jll[j], env->memb_idx_ju[j], env->memb_idx_jb[j], env->memb_idx_jbb[j]);
		new_env->cell_y[j] =celly+ CONST::dt_cell*bend_force_sqr(ny, my, ipn, ipm, inv_dn, inv_dm, j, env->memb_idx_jr[j], env->memb_idx_jl[j], env->memb_idx_jll[j], env->memb_idx_ju[j], env->memb_idx_jb[j], env->memb_idx_jbb[j]);
		new_env->cell_z[j] =cellz+ CONST::dt_cell*bend_force_sqr(nz, mz, ipn, ipm, inv_dn, inv_dm, j, env->memb_idx_jr[j], env->memb_idx_jl[j], env->memb_idx_jll[j], env->memb_idx_ju[j], env->memb_idx_jb[j], env->memb_idx_jbb[j]);
	}
	
}

void interac_wall(double old_zz, double r, double *zz)
{
	double inv_zz = 1 / old_zz;
	double rinvzz = inv_zz*r;
	double rinvzzSq = rinvzz*rinvzz;
	double dum_lj6 = rinvzzSq*rinvzzSq*rinvzzSq;
	*zz +=2.0* CONST::dt_ca*CONST::eps_m * dum_lj6 * (dum_lj6 - 1.0) *inv_zz;
}



void memb_memb_interact(const ENV* env, ENV* new_env) {
	for (int j = 0; j < env->NMEMB; j++) {
		double cellx, celly, cellz, diffx, diffy, diffz, distlj_inv, rad_sum, cr_dist, cr_dist_inv, lambda_dist, cr_lj, dum_lj6, dum_ljmain;
		double cell_con_x, cell_con_y, cell_con_z;
		double dumxx, dumyy, dumzz;
		cellx = env->cell_x[j];
		celly = env->cell_y[j];
		cellz = env->cell_y[j];
		//半径の方がz位置(高さ)より大きい->z=0と接している
		if (env->cell_z[j] < env->cell_radius[j]) {
			interac_wall(env->cell_z[j], env->cell_radius[j], &(new_env->cell_z[j]));
		}

		for (int k = 0; k < env->cell_connected_num[j]; k++) {
			int connected_index = env->cell_connected_index[j][k];
			if (env->cell_state[connected_index] != MEMB || j>connected_index)continue;
			if (k >= CONST::N_NEIGHBOR) {
				printf("error: membrane particles should be labeled k < NN\n");
				exit(1);
			}
			cell_con_x = env->cell_x[connected_index];
			cell_con_y = env->cell_y[connected_index];
			cell_con_z = env->cell_z[connected_index];
			diffx = cellx-ret_x(cellx, cell_con_x);
			diffy = celly-ret_y(celly , cell_con_y);
			diffz = cellz - cell_con_z;

			distlj_inv = t_sqrtD(diffx*diffx + diffy*diffy + diffz*diffz);
			rad_sum = (env->cell_radius[j] + env->cell_radius[connected_index]);
			cr_dist = rad_sum*CONST::P_MEMB;
			cr_dist_inv = 1 / cr_dist;
			 lambda_dist = (1.0 + CONST::P_MEMB)*rad_sum;
			 cr_lj = cr_dist*distlj_inv;
			 dum_lj6,dum_ljmain;
			if (cr_lj > 1) { //distlj < cr_dist
				double cr_ljSq = cr_lj*cr_lj;
				dum_lj6 = cr_ljSq*cr_ljSq*cr_ljSq;
				dum_ljmain = 4 * CONST::eps_m*dum_lj6*(dum_lj6 - 1.0)*distlj_inv*distlj_inv;
			}
			else if (rad_sum*distlj_inv > 1) { //distlj<ra+rb
				dum_ljmain = CONST::DER_DER_CONST*(distlj_inv - cr_dist_inv);
			}
			else {
				double lamb_distlj_inv = 1 / (lambda_dist*distlj_inv - 1.0);
				double cr_lj_lamb = cr_lj *lamb_distlj_inv;
				double cr_lj_lambSq = cr_lj_lamb*cr_lj_lamb;
				dum_lj6 = cr_lj_lambSq*cr_lj_lambSq*cr_lj_lambSq;

				dum_ljmain = CONST::DER_DER_CONST*(CONST::P_MEMB - 1.0)*cr_dist_inv
					- 4.0*CONST::eps_m*dum_lj6*(dum_lj6 - 1.0)*distlj_inv*distlj_inv*lamb_distlj_inv;

			}

			if (dum_ljmain*dum_ljmain > 10000) {
				printf("error: interaction der-der too strong\n");
				//printf("id=%d state=%d l=%d state=%d xx=%f yy=%f zz=%f, dum_ljmain=%f\n", j, old_state[j], l, old_state[l], xx[j], yy[j], zz[j], dum_ljmain);
				exit(1);
			}
			dumxx = CONST::dt_cell*dum_ljmain*diffx; //calced in dist3 inv
			dumyy = CONST::dt_cell*dum_ljmain*diffy;
			dumzz = CONST::dt_cell*dum_ljmain*diffz;
			new_env->cell_x[j] += dumxx;
			new_env->cell_y[j] += dumyy;
			new_env->cell_z[j] += dumzz;
			new_env->cell_x[connected_index] -= dumxx;
			new_env->cell_y[connected_index] -= dumyy;
			new_env->cell_z[connected_index] -= dumzz;
		}
	}
}

void interact_lj(int j, int connected_index, const ENV* env, ENV* new_env) {
	double cellx = env->cell_x[j];
	double celly = env->cell_y[j];
	double cellz = env->cell_z[j];
	double diffx = cellx-ret_x(cellx, env->cell_x[connected_index]);
	double diffy = celly-ret_y(celly, env->cell_y[connected_index]);
	double diffz = cellz - env->cell_z[connected_index];
	double distlj_inv = t_sqrtD(diffx*diffx + diffy*diffy + diffz*diffz);
	double LJ1 = (env->cell_radius[j] + env->cell_radius[connected_index])*distlj_inv;
	double LJ1Sq = LJ1*LJ1;
	double dum_lj6 = LJ1Sq*LJ1Sq*LJ1Sq;
	double dum_ljmain = 4.0*CONST::eps_m*dum_lj6*(dum_lj6 - 1.0)*distlj_inv*distlj_inv;

	double dumxx = CONST::dt_cell*dum_ljmain*diffx;
	double dumyy= CONST::dt_cell*dum_ljmain*diffy;
	double dumzz = CONST::dt_cell*dum_ljmain*diffz;
	new_env->cell_x[j] += dumxx;
	new_env->cell_y[j] += dumyy;
	new_env->cell_z[j] += dumzz;

	new_env->cell_x[connected_index] -= dumxx;
	new_env->cell_y[connected_index] -= dumyy;
	new_env->cell_z[connected_index] -= dumzz;
}
void der_der_interact_DER(int cell_idx,int cell_con_idx,double cellx, double celly, double cellz,
	double cell_con_x, double cell_con_y, double cell_con_z,
	double cell_rad, double cell_con_rad, ENV* new_env) {
	double dum_lj6, dum_ljmain;
	double diffx = cellx-ret_x(cellx, cell_con_x);
	double diffy = celly-ret_y(celly, cell_con_y);
	double diffz = cellz - cell_con_z;
	double distlj_inv = t_sqrtD(diffx*diffx + diffy*diffy + diffz*diffz);
	double distlj_delta_inv = distlj_inv / (1.0 + CONST::delta_R*distlj_inv);
	double rad_sum = cell_rad + cell_con_rad;
	double LJ1 = rad_sum*distlj_delta_inv;
	if (rad_sum*distlj_inv > 1) {
		dum_ljmain = CONST::ljp2;
	}
	else {
		double LJ1Sq = LJ1*LJ1;
		dum_lj6 = LJ1Sq*LJ1Sq*LJ1Sq;
		dum_ljmain = 4.0*CONST::eps_m*dum_lj6*(dum_lj6 - 1.0)*(LJ1>1 ? distlj_delta_inv : distlj_inv)*distlj_inv + CONST::ljp2;
	}
	double dumxx = CONST::dt_cell*dum_ljmain*diffx;
		double dumyy = CONST::dt_cell*dum_ljmain*diffy;
	double dumzz = CONST::dt_cell*dum_ljmain*diffz;
	new_env->cell_x[cell_idx] += dumxx;
	new_env->cell_y[cell_idx] += dumyy;
	new_env->cell_z[cell_idx] += dumzz;

	new_env->cell_x[cell_con_idx] -= dumxx;
	new_env->cell_y[cell_con_idx] -= dumyy;
	new_env->cell_z[cell_con_idx] -= dumzz;
}
void der_der_interact(const ENV* env, ENV* new_env) {
	for (int j = env->NMEMB; j < env->NMEMB + env->NDER; j++) {
		//半径の方がz位置(高さ)より大きい->z=0と接している
		if (env->cell_z[j] < env->cell_radius[j]) {
			interac_wall(env->cell_z[j], env->cell_radius[j], &(new_env->cell_z[j]));
		}

		for (int k = 0; k < env->cell_connected_num[j]; k++) {
			int connected_index = env->cell_connected_index[j][k];
			CELL_STATE state = env->cell_state[connected_index];

			assert(state&(ALIVE | DEAD | AIR | DER | MEMB | FIX | MUSUME));

			if (state&(MEMB | FIX | MUSUME | ALIVE | DEAD | AIR)) {
				interact_lj(j, connected_index, env, new_env);
			}
			else if (state == DER) {
				der_der_interact_DER(j, connected_index, env->cell_x[j], env->cell_y[j], env->cell_z[j],
					env->cell_x[connected_index], env->cell_y[connected_index], env->cell_z[connected_index],
					env->cell_radius[j], env->cell_radius[connected_index], new_env);
			}
		}
	}
}

inline double adhesion(double distlj, double r_sum, double spring_const)
{

	double LJ2_m1 = distlj / r_sum-1;
	return (LJ2_m1+1 > CONST::LJ_THRESH ?
		0.0 :
		-(spring_const / distlj)*LJ2_m1
		* (1 - LJ2_m1*LJ2_m1 / ((CONST::LJ_THRESH - 1.0)*(CONST::LJ_THRESH - 1.0))));
}

void other_interact_ALIVE_DEAD_AIR(int cell_idx, int cell_con_idx, double cellx, double celly, double cellz,
	double cell_con_x, double cell_con_y, double cell_con_z,
	double cell_rad, double cell_con_rad,double cell_agek,double cell_con_agek, ENV* new_env) {
	double dum_lj6, dum_ljmain;
	double LJ1, LJ1Sq;
	double diffx = cellx-ret_x(cellx, cell_con_x);
	double diffy = celly-ret_y(celly, cell_con_y);
	double diffz = cellz - cell_con_z;
	double distlj_inv = t_sqrtD(diffx*diffx + diffy*diffy + diffz*diffz);
	double rad_sum = cell_rad + cell_con_rad;
	
	double spring_force;
	if (rad_sum*distlj_inv > 1) {
		LJ1 = distlj_inv*rad_sum;
		LJ1Sq = LJ1*LJ1;
		dum_lj6 = LJ1Sq*LJ1Sq*LJ1Sq;
		dum_ljmain = 4.0*CONST::eps_m*dum_lj6*(dum_lj6 - 1.0)*distlj_inv*distlj_inv;
	}
	else {
		//attr
		spring_force = cell_agek > CONST::THRESH_SP && cell_con_agek > CONST::THRESH_SP ?
			CONST::K_TOTAL :
			CONST::K_DESMOSOME;
		//adhe
		//check precision...
		dum_ljmain = adhesion(1 / distlj_inv, rad_sum, CONST::K_DESMOSOME);
	}
	double dumxx = CONST::dt_cell*dum_ljmain*diffx;
	double dumyy = CONST::dt_cell*dum_ljmain*diffy;
	double dumzz = CONST::dt_cell*dum_ljmain*diffz;
	new_env->cell_x[cell_idx] += dumxx;
	new_env->cell_y[cell_idx] += dumyy;
	new_env->cell_z[cell_idx] += dumzz;

	new_env->cell_x[cell_con_idx] -= dumxx;
	new_env->cell_y[cell_con_idx] -= dumyy;
	new_env->cell_z[cell_con_idx] -= dumzz;
}

int find_dermis_index(const ENV* env, int cell_idx,double cellx,double celly,double cellz) {
	double dist1 = (CONST::Lx*CONST::Lx);
	int dermis_idx = -1;
	int cell_con_num = env->cell_connected_num[cell_idx];
	double diffx, diffy, diffz,distSq;
	for (int k = 0; k < cell_con_num; k++) {
		int connected_index = env->cell_connected_index[cell_idx][k];
		diffx = cellx-ret_x(cellx, env->cell_x[connected_index]);
		diffy = celly-ret_y(celly , env->cell_y[connected_index]);
		diffz = cellz - env->cell_z[connected_index];
		distSq = diffx*diffx + diffy*diffy + diffz*diffz;
		if (env->cell_state[connected_index] == MEMB && distSq < dist1) {
			dist1 = distSq;
			dermis_idx = connected_index;
		}

	}
	return dermis_idx;
}

void other_interact_stem_memb(int cell_idx, int cell_con_idx, double cellx, double celly, double cellz,
	double cell_con_x, double cell_con_y, double cell_con_z,
	double cell_rad, double cell_con_rad, double spring_force, CELL_STATE cell_state, ENV* new_env) {
	double dum_lj6, dum_ljmain;
	double LJ1, LJ1Sq;
	double diffx = cellx-ret_x(cellx, cell_con_x);
	double diffy = celly-ret_y(celly, cell_con_y);
	double diffz = cellz - cell_con_z;
	double distlj_inv = t_sqrtD(diffx*diffx + diffy*diffy + diffz*diffz);
	double rad_sum = cell_rad + cell_con_rad;

	if (rad_sum*distlj_inv > 1) {
		LJ1 = distlj_inv*rad_sum;
		LJ1Sq = LJ1*LJ1;
		dum_lj6 = LJ1Sq*LJ1Sq*LJ1Sq;
		dum_ljmain = 4.0*CONST::eps_m*dum_lj6*(dum_lj6 - 1.0)*distlj_inv*distlj_inv;
	}
	else {
		dum_ljmain = cell_state == FIX ? -(spring_force)*(1 / rad_sum - distlj_inv) : adhesion(1 / distlj_inv, rad_sum, spring_force);

	}
	assert(dum_ljmain*dum_ljmain <= 10000);
	double dumxx = CONST::dt_cell*dum_ljmain*diffx;
	double dumyy = CONST::dt_cell*dum_ljmain*diffy;
	double dumzz = CONST::dt_cell*dum_ljmain*diffz;
	new_env->cell_x[cell_idx] += dumxx;
	new_env->cell_y[cell_idx] += dumyy;
	new_env->cell_z[cell_idx] += dumzz;

	new_env->cell_x[cell_con_idx] -= dumxx;
	new_env->cell_y[cell_con_idx] -= dumyy;
	new_env->cell_z[cell_con_idx] -= dumzz;
}

void other_interact_stem_stem(int cell_idx, int cell_con_idx, double cellx, double celly, double cellz,
	double cell_con_x, double cell_con_y, double cell_con_z,
	double cell_rad, double cell_con_rad,int pair_index, CELL_STATE cell_state, ENV* new_env) {
	double dum_lj6, dum_ljmain;
	double diffx = cellx - ret_x(cellx, cell_con_x);
	double diffy = celly - ret_y(celly, cell_con_y);
	double diffz = cellz - cell_con_z;
	if (pair_index == cell_con_idx) {
		dum_ljmain = 0;
	}
	else {
		double LJ1, LJ1Sq;
		
		
		double distlj_inv = t_sqrtD(diffx*diffx + diffy*diffy + diffz*diffz);
		double rad_sum = cell_rad + cell_con_rad;

		if (rad_sum*distlj_inv > 1) {
			LJ1 = distlj_inv*rad_sum;
			LJ1Sq = LJ1*LJ1;
			dum_lj6 = LJ1Sq*LJ1Sq*LJ1Sq;
			dum_ljmain = 4.0*CONST::eps_m*dum_lj6*(dum_lj6 - 1.0)*distlj_inv*distlj_inv;
		}
		else {
			dum_ljmain = adhesion(1 / distlj_inv, rad_sum, CONST::K_DESMOSOME);

		}
	}

	assert(dum_ljmain*dum_ljmain <= 10000);

	double dumxx = CONST::dt_cell*dum_ljmain*diffx;
	double dumyy = CONST::dt_cell*dum_ljmain*diffy;
	double dumzz = CONST::dt_cell*dum_ljmain*diffz;
	new_env->cell_x[cell_idx] += dumxx;
	new_env->cell_y[cell_idx] += dumyy;
	new_env->cell_z[cell_idx] += dumzz;

	new_env->cell_x[cell_con_idx] -= dumxx;
	new_env->cell_y[cell_con_idx] -= dumyy;
	new_env->cell_z[cell_con_idx] -= dumzz;
}
void other_interact(const ENV* env, ENV* new_env) {
//	double diffx, diffy, diffz, dumxx, dumyy, dumzz,distlj_inv,rad_sum,LJ1,LJ1Sq,dum_lj6,dum_ljmain;
	double cellx, celly, cellz,  cell_agek,cell_rad;
	int cell_pair_idx,cell_div_times,cell_con_num,dermis_index;
	double spring_force;
	CELL_STATE cell_state,cell_pair_state;
	for (int j = env->NMEMB + env->NDER; j < env->current_cell_num; j++) {

		cellx = env->cell_x[j];
		celly = env->cell_y[j];
		cellz = env->cell_z[j];
		cell_agek = env->cell_agek[j];
		cell_rad = env->cell_radius[j];
		cell_state = env->cell_state[j];
		cell_pair_idx = env->cell_div_pair_index[j];
		cell_div_times = env->cell_div_times[j];
		cell_pair_state = cell_pair_idx == -1 ? CELL_UNDEF : env->cell_state[cell_pair_idx];
		cell_con_num = env->cell_connected_num[j];

		switch (env->cell_state[j]) {
		case DISA:
			break;

		case ALIVE:case DEAD:case AIR:
			for (int k = 0; k < cell_con_num; k++) {
				int connected_index = env->cell_connected_index[j][k];
				CELL_STATE state = env->cell_state[connected_index];
				if (state==DER || ((state &(ALIVE|DEAD|AIR)) && connected_index>j)) continue;

				other_interact_ALIVE_DEAD_AIR(j, connected_index, env->cell_x[j], env->cell_y[j], env->cell_z[j],
					env->cell_x[connected_index], env->cell_y[connected_index], env->cell_z[connected_index],
					env->cell_radius[j], env->cell_radius[connected_index],cell_agek, env->cell_agek[connected_index], new_env);

			}
			break;

		case FIX:case MUSUME:
			spring_force = 0;
			if (cell_state == FIX || (cell_state == MUSUME && cell_pair_state == FIX)) {
				spring_force = CONST::Kspring;
			}
			else if (cell_div_times > 0) {
				spring_force = CONST::Kspring_d;
			}
			dermis_index = find_dermis_index(env, j, cellx, celly, cellz);

			for (int k = 0; k < cell_con_num; k++) {
				int connected_index = env->cell_connected_index[j][k];
				CELL_STATE cell_con_state = env->cell_state[connected_index];
				assert(cell_con_state&(ALIVE | DEAD | AIR | DER | MEMB | FIX | MUSUME));
				if (cell_con_state == MEMB) {
					if (connected_index == dermis_index) {
						other_interact_stem_memb(j, connected_index,
							cellx, celly, cellz,
							env->cell_x[connected_index], env->cell_y[connected_index], env->cell_z[connected_index],
							cell_rad, env->cell_radius[connected_index],
							spring_force, cell_state, new_env);
					}
					else {
						interact_lj(j, connected_index, env, new_env);
					}
				}
				else if (cell_con_state &(FIX | MUSUME)) {
					if (j < connected_index) {
						other_interact_stem_stem(j, connected_index,
							cellx, celly, cellz,
							env->cell_x[connected_index], env->cell_y[connected_index], env->cell_z[connected_index],
							cell_rad, env->cell_radius[connected_index], 
							cell_pair_idx, cell_state, new_env);
					}
				}

				
			}
			break; //coffee
		default:
			break;
		}
	}


}

void pair_interact(const ENV* env, ENV* new_env) {
	int pair_count = env->pair_count;
	int less_idx, large_idx;
	double diffx, diffy, diffz,force,dumxx,dumyy,dumzz,dist_inv, nat_spr_len;
	double le_x, le_y, le_z, la_x, la_y, la_z;
	for (int k = 0; k < pair_count; k++) {
		less_idx = env->cell_pair_lesser[k];
		large_idx = env->cell_pair_larger[k];
		nat_spr_len = env->cell_spring_len[less_idx];

		le_x = env->cell_x[less_idx]; le_y = env->cell_y[less_idx]; le_z = env->cell_z[less_idx];
		la_x = env->cell_x[large_idx]; la_y = env->cell_y[large_idx]; la_z = env->cell_z[large_idx];

		diffx = le_x-ret_x(le_x, la_x); diffy = le_y-ret_y(le_y, la_x); diffz = le_z - la_z;

		dist_inv = t_sqrtD(diffx*diffx + diffy*diffy + diffz*diffz);

		//cjeck precision
		force = CONST::Kspring_division*(1 / dist_inv - nat_spr_len);

		dumxx = CONST::dt_cell*force*diffx*dist_inv;
		dumyy = CONST::dt_cell*force*diffy*dist_inv;
		dumzz = CONST::dt_cell*force*diffz*dist_inv;
		new_env->cell_x[less_idx] += dumxx;
		new_env->cell_y[less_idx] += dumyy;
		new_env->cell_z[less_idx] += dumzz;

		new_env->cell_x[large_idx] -= dumxx;
		new_env->cell_y[large_idx] -= dumyy;
		new_env->cell_z[large_idx] -= dumzz;
	}
}

void fixed_memb_pos(const ENV* env, ENV* new_env) {
	if (CONST::KEY_DERMAL_CHANGE == 0) {
		memcpy(new_env->cell_x, env->cell_x, CONST::MAX_CELL_NUM*sizeof(double));
		memcpy(new_env->cell_y, env->cell_y, CONST::MAX_CELL_NUM*sizeof(double));
		memcpy(new_env->cell_z, env->cell_z, CONST::MAX_CELL_NUM*sizeof(double));
	}
}
//do not execute this in parallel
void calc_dermis_normal(double cellx, double celly, double cellz, double der_x, double der_y, double der_z,
	double *nx, double *ny, double *nz, double *mx, double *my, double *mz) {

	//original code does something wasteful

	*mx = ret_x(cellx, der_x);
	*my = ret_y(celly, der_y);
	*mz = der_z;
	double nx1 = cellx - *mx;
	double ny1 = celly - *my;
	double nz1 = cellz - *mz;

	double norm_inv = t_sqrtD(nx1*nx1 + ny1*ny1 + nz1*nz1);
	*nx = nx1*norm_inv;
	*ny = ny1*norm_inv;
	*nz = nz1*norm_inv;

}

void div_direction(double cellx,double celly,double cellz,
	double der_x,double der_y,double der_z,
	double *nx, double *ny, double *nz) {
	double lx, ly, lz; double mx, my, mz;
	double nx1, ny1, nz1;
	double rand_theta, rand_phi, sr1;
	double randnx, randny, randnz, sum;
	calc_dermis_normal(cellx, celly, cellz, der_x, der_y, der_z, &lx, &ly, &lz, &mx, &my, &mz);
	
	do {
		rand_theta = M_PI*genrand_res53();
		rand_phi = 2 * M_PI*genrand_res53();
		sr1 = sin(rand_theta);
		randnx = sr1*cos(rand_phi);
		randny = sr1*sin(rand_phi);
		randnz = cos(rand_theta);

		nx1 = ly*randnz - lz*randny;
		ny1 = lz*randnx - lx*randnz;
		nz1 = lx*randny - ly*randnx;
	} while ((sum = t_sqrtD(nx1*nx1 + ny1*ny1 + nz1*nz1)) > 1.0e+7);

	*nx = nx1*sum;
	*ny = ny1*sum;
	*nz = nz1*sum;
}

void set_new_state(CELL_STATE* old, CELL_STATE* new_state) {
	*old = *new_state = MUSUME;
}

inline double ageb_const(CELL_STATE state, double ca, int orig_stem_idx) {
	double subu = ca - CONST::u0;
	return (orig_stem_idx < CONST::MALIGNANT_NUM ? CONST::accel_div : 1)*CONST::eps_kb*(CONST::u0 + CONST::alpha_b*(subu>0 ? subu : 0));
}

void cell_div(ENV* env,ENV *update_env,int cell_idx,int given_divtimes,int dermis_index){
	//書き換えの順番は要チェック

	//分裂中
	assert(dermis_index != -1);
	if (env->cell_div_pair_index[cell_idx] != -1)return;
	CELL_STATE cell_state = env->cell_state[cell_idx];
	double par_inv,div_gamma;
	int new_cell_idx;
	par_inv = cell_state == MUSUME ?
		CONST::dt_cell*CONST::eps_kb*CONST::u0*CONST::stoch_div_time_ratio_inv*CONST::agki_max_inv :
		CONST::dt_cell*CONST::eps_kb*CONST::u0*CONST::stoch_div_time_ratio_inv*CONST::agki_max_fix_inv; //constexpr

	div_gamma = env->origin_stem_cell_index[cell_idx] < CONST::MALIGNANT_NUM ?
		CONST::accel_div*par_inv : par_inv;

	if (CONST::STOCHASTIC && genrand_res53() > div_gamma)return;

	update_env->born_num = env->born_num + 1;
	
	if (env->vacant_num == 0) {
		new_cell_idx = env->current_cell_num;

		update_env->current_cell_num = env->current_cell_num + 1;
		//errorcheck
	}
	else {
		new_cell_idx = env->vanished_cell_index[env->vacant_num];
		update_env->vacant_num = env->vacant_num + 1;
	}
	env->cell_generated[new_cell_idx] = true; //forbid overwriting


	update_env->cell_div_pair_index[cell_idx] = new_cell_idx;
	update_env->cell_div_pair_index[new_cell_idx] = cell_idx;

	assert(new_cell_idx != cell_idx);

	update_env->cell_pair_lesser[update_env->next_pair_index] = min(cell_idx, new_cell_idx);
	update_env->cell_pair_larger[update_env->next_pair_index] = max(cell_idx, new_cell_idx);
	update_env->cell_spring_len[cell_idx] = update_env->cell_spring_len[new_cell_idx] = CONST::delta_L;
	update_env->next_pair_index++;

	double x_orig = env->cell_x[cell_idx];
	double y_orig = env->cell_y[cell_idx];
	double z_orig = env->cell_z[cell_idx];

	double rh = CONST::delta_L*0.5;

	update_env->cell_radius[new_cell_idx] = env->cell_radius[cell_idx];
	update_env->cell_div_times[new_cell_idx] = given_divtimes;

	if (cell_state == MUSUME) {
		update_env->cell_div_times[new_cell_idx]--;
		update_env->cell_div_times[cell_idx] = env->cell_div_times[cell_idx] - 1;
	}
	double dir_nx, dir_ny, dir_nz;
	div_direction(env->cell_x[cell_idx], env->cell_y[cell_idx], env->cell_z[cell_idx],
		env->cell_x[dermis_index], env->cell_y[dermis_index], env->cell_z[dermis_index],
		&dir_nx, &dir_ny, &dir_nz);
	dir_nx *= rh; dir_ny *= rh; dir_nz *= rh;
	update_env->cell_x[cell_idx] = x_orig + dir_nx;
	update_env->cell_y[cell_idx] =y_orig + dir_ny;
	update_env->cell_z[cell_idx] = z_orig + dir_nz;

	update_env->cell_x[new_cell_idx] = x_orig - dir_nx;
	update_env->cell_y[new_cell_idx] = y_orig - dir_ny;
	update_env->cell_z[new_cell_idx] = z_orig - dir_nz;

	//set_new_state(update_env->cell_state[)
	update_env->cell_state[new_cell_idx] = MUSUME;

	update_env->cell_c[new_cell_idx] = update_env->cell_c[cell_idx];
	update_env->cell_ca_ave[new_cell_idx] = update_env->cell_ca_ave[cell_idx];
	update_env->cell_h[new_cell_idx] = update_env->cell_h[cell_idx];
	update_env->cell_P[new_cell_idx] = update_env->cell_P[cell_idx];
	update_env->origin_stem_cell_index[new_cell_idx] = update_env->origin_stem_cell_index[cell_idx];

	update_env->cell_ageb[new_cell_idx] = 0.0;
	update_env->cell_agek[new_cell_idx] = 0.0;
	update_env->cell_fat[new_cell_idx] = 0.0;
	update_env->cell_Vc[new_cell_idx] = 0.0;


    //double par=env->cell_state[target_idx]==MUSUME?CONST::ag
}
inline double k_lipid_release(double u, double age) {

	return 0.25*CONST::lipid_rel*(1 + tanh(100 * (u - CONST::ubar)))*(1 + tanh((age - CONST::THRESH_SP)*CONST::delta_lipid_inv));
}
void state_renew(ENV* env,ENV* new_env) {
	int dermis_index;
	double cellx, celly, cellz;
	CELL_STATE cst;
	for (int j = env->current_cell_num - 1; j >= env->NMEMB + env->NDER; j--) {
		cellx = env->cell_x[j];
		celly = env->cell_y[j];
		cellz = env->cell_z[j];
		if (env->cell_generated[j])continue; //forbid writing
		
		switch (cst=env->cell_state[j])
		{
		case FIX:
			dermis_index = find_dermis_index(env, j, cellx, celly, cellz);
			assert(dermis_index != -1);
			if (env->cell_ageb[j] >= CONST::agki_max_fix*(1.0 - CONST::stoch_div_time_ratio)) {
				cell_div(env, new_env, j, CONST::div_max, dermis_index);
			}
			else {
				new_env->cell_ageb[j] =
					env->cell_ageb[j] + CONST::dt_cell*ageb_const(cst, env->cell_ca_ave[j], env->origin_stem_cell_index[j]);
			}
			break;
		case MUSUME:
			dermis_index = find_dermis_index(env, j, cellx, celly, cellz);
			if (dermis_index == -1 && env->cell_div_pair_index[j] == -1) {
				if (c_SYSTEM == WHOLE) {
					new_env->cell_state[j] = ALIVE;

				}
				else if(c_SYSTEM==BASAL) {
					new_env->cell_state[j] = DISA;
					new_env->vacant_num = env->vacant_num + 1;
					new_env->vanished_cell_index[new_env->vacant_num] = j;
					new_env->disap_num = env->disap_num + 1;
				}
			}
			else {
				if (env->cell_div_times[j]>0&&env->cell_ageb[j] >= CONST::agki_max_fix*(1.0 - CONST::stoch_div_time_ratio)) {
					cell_div(env, new_env, j, CONST::div_max, dermis_index);
				}
				else {
					new_env->cell_ageb[j] =
						env->cell_ageb[j] + CONST::dt_cell*ageb_const(cst, env->cell_ca_ave[j], env->origin_stem_cell_index[j]);
				}
			}
			break;

		case DEAD:case AIR:
			if (env->cell_agek[j] >= CONST::ADHE_CONST && env->cell_connected_num[j] <= CONST::Nc) {
				new_env->cell_state[j] = DISA;
				//del
				new_env->vacant_num = env->vacant_num + 1;
				new_env->vanished_cell_index[new_env->vacant_num] = j;
				new_env->disap_num = env->disap_num + 1;
			}
			else {
				new_env->cell_agek[j] = env->cell_agek[j] +
					CONST::dt_cell*ageb_const(cst, CONST::u0, env->origin_stem_cell_index[j]);
			}
			break;
		case ALIVE:
			if (env->cell_agek[j] >= CONST::THRESH_DEAD) {
				new_env->cell_state[j] = DEAD;
				new_env->sw = env->sw + 1;
			}
			else {
				new_env->cell_agek[j] =
					env->cell_agek[j] + CONST::dt_cell*ageb_const(cst, env->cell_ca_ave[j], env->origin_stem_cell_index[j]);
				double klrel = k_lipid_release(env->cell_ca_ave[j], env->cell_agek[j]);
				double tmp = klrel*env->cell_Vc[j];
				new_env->cell_Vc[j] = env->cell_Vc[j] + CONST::dt_cell*klrel*(1.0 - 2.0*env->cell_Vc[j]);
				new_env->cell_fat[j] = env->cell_fat[j] + CONST::dt_cell*tmp;
			}
			break;
		default:
			break;
		}
	}
}
void moveover(const ENV* env,ENV* new_env, int nunpaired)
{
	int j, jj; 
	std::vector<int> tmp(env->pair_count,0),tmp2(env->pair_count,0);
	for (j = 0, jj = 0; j < env->pair_count; j++) {
		if (env->cell_pair_lesser[j] != -1) {
			new_env->cell_pair_lesser[jj] = env->cell_pair_lesser[j];
			new_env->cell_pair_larger[jj] = env->cell_pair_larger[j];
			jj++;
		}
	}
	for (j = env->pair_count - nunpaired; j < env->pair_count; j++) {
		new_env->cell_pair_lesser[j] = -1;
		new_env->cell_pair_larger[j] = -1;
	}

}

void pair_disperse(const ENV* env, ENV* new_env) {
	int unpaired_count = 0;
	for (int k = 0; k < env->pair_count; k++) {
		int less_idx = env->cell_pair_lesser[k];
		int large_idx = env->cell_pair_larger[k];

		if (env->cell_spring_len[less_idx] < 2.0*env->cell_radius[less_idx]) {
			new_env->cell_spring_len[large_idx]=new_env->cell_spring_len[less_idx] =
				env->cell_spring_len[less_idx] + CONST::dt_ca*CONST::eps_L;
		}
		else {
			double diffx = env->cell_x[less_idx] - ret_x(env->cell_x[less_idx],env->cell_x[large_idx]);
			double diffy = env->cell_y[less_idx] - ret_y(env->cell_y[less_idx],env->cell_y[large_idx]);
			double diffz = env->cell_z[less_idx] - env->cell_z[large_idx];
			double rad_sum = env->cell_radius[less_idx] + env->cell_radius[large_idx];
			if (diffx*diffx + diffy*diffy + diffz*diffz>CONST::div_fin_coef*CONST::div_fin_coef*rad_sum*rad_sum) {
				new_env->cell_spring_len[large_idx] = new_env->cell_spring_len[less_idx] = 0;
				new_env->cell_div_pair_index[less_idx] = -1;
				new_env->cell_div_pair_index[large_idx] = -1;

				new_env->cell_pair_lesser[k] = -1;
				new_env->cell_pair_larger[k] = -1;
				unpaired_count++;
			}
		}
	}
	if (unpaired_count > 0) {
		moveover(env, new_env, unpaired_count);
		new_env->pair_count = env->pair_count - unpaired_count;
	}
}

void pos_fix_periodic_cond(ENV* env)
{
	int j;

	for (j = 0; j < env->pair_count; j++) {
		if (env->cell_x[j] > CONST::Lx) {
			env->cell_x[j] -= CONST::Lx;
		}
		else if (env->cell_x[j] < 0.0) {
			env->cell_x[j] += + CONST::Lx;
		}

		if (env->cell_y[j] > CONST::Ly) {
			env->cell_y[j] -= CONST::Ly;
		}
		else if (env->cell_y[j] < 0.0) {
			env->cell_y[j] += +CONST::Ly;
		}
	}
}

void connect_area_calc(const ENV* env,
	int area_cell_idx[CONST::ANX][CONST::ANY][CONST::ANZ][CONST::cell_on_lat_max],
	int area_cell_num[CONST::ANX][CONST::ANY][CONST::ANZ]) {
	std::fill(&(area_cell_num[0][0][0]), &(area_cell_num[0][0][0])+CONST::ANX*CONST::ANY*CONST::ANZ, 0);
	int cnum = env->current_cell_num;
	int aix, aiy, aiz;
	for (int j = 0; j < cnum; j++) {
		if (env->cell_state[j] != DISA) {
			aix = (int)(ret_x(CONST::Lx *0.5, env->cell_x[j]) * CONST::AREA_GRID_inv);
			aiy = (int)(ret_y(CONST::Ly *0.5, env->cell_y[j]) * CONST::AREA_GRID_inv);
			aiz = (int)((env->cell_z[j] < 0. ? 0 : env->cell_z[j]) * CONST::AREA_GRID_inv);

			assert(!(aix >= CONST::ANX || aiy >= CONST::ANY || aiz >= CONST::ANZ || aix < 0 || aiy < 0 || aiz < 0));

			area_cell_idx[aix][aiy][aiz][area_cell_num[aix][aiy][aiz]] = j;
			area_cell_num[aix][aiy][aiz]++;
#ifndef TEST
			assert(area_cell_num[aix][aiy][aiz] < CONST::cell_on_lat_max);
#endif // TEST

			
		}
	}
}

void connect_lj(ENV* env) {
	int jj, kk,jl,jr,jb,ju,ii;
	int anx, any, anz;
	int aix, aiy, aiz;
	static int indx_bak[CONST::MAX_CELL_NUM];
	static int area_cell_idx[CONST::ANX][CONST::ANY][CONST::ANZ][CONST::cell_on_lat_max];
	static int area_cell_num[CONST::ANX][CONST::ANY][CONST::ANZ];
	CELL_STATE state;
	connect_area_calc(env, area_cell_idx, area_cell_num);
	//optimize point
	for (int j = 0; j < env->NMEMB; j++) {
		//NMEMB 横につなぐ？
		env->cell_connected_num[j] = 0; 
		jj = j%CONST::NMX;
		kk = j / CONST::NMX;
		if (jj == 0) {
			jl = j + CONST::NMX - 1;
		}
		else {
			jl = j - 1;
		}

		if (jj == CONST::NMX - 1) {
			jr = j - (CONST::NMX - 1);
		}
		else {
			jr = j + 1;
		}

		if (kk == 0) {
			jb = j + (CONST::NMX*CONST::NMY - CONST::NMX);

		}
		else {
			jb = j - CONST::NMX;
		}

		if (kk == CONST::NMY - 1) {
			ju = j - (CONST::NMX*CONST::NMY - CONST::NMX);
		}
		else {
			ju = j + CONST::NMX;
		}

		env->cell_connected_index[j][env->cell_connected_num[j]] = jl; env->cell_connected_num[j]++;
		env->cell_connected_index[j][env->cell_connected_num[j]] = jr; env->cell_connected_num[j]++;
		env->cell_connected_index[j][env->cell_connected_num[j]] = jb; env->cell_connected_num[j]++;
		env->cell_connected_index[j][env->cell_connected_num[j]] = ju; env->cell_connected_num[j]++;

		indx_bak[j] = env->cell_connected_num[j];
	}
	double cellx, celly, cellz;
	bool membflg = false;
#pragma omp parallel for
	for (int i = 0; i < env->current_cell_num; i++) {
		state = env->cell_state[i];
		if (state == DISA) continue;
		membflg = state == MEMB;
		cellx = env->cell_x[i];
		celly = env->cell_y[i];
		cellz = env->cell_z[i];
		env->cell_connected_num[i] = (state == MEMB ? indx_bak[i] : 0);
		anx = (int)(cellx *CONST::AREA_GRID_inv);
		any = (int)(celly * CONST::AREA_GRID_inv);
		anz = (int)(cellz * CONST::AREA_GRID_inv);

		assert(!(anx >= CONST::ANX || any >= CONST::ANY || anz >= CONST::ANZ || anx < 0 || any < 0 || anz < 0));
	
		ii = 2;
		int o;
		double cell_rad = env->cell_radius[i];
		double diffx, diffy, diffz;
		for (int j = anx - ii; j <= anx + ii; j++) {
			aix = (j + CONST::ANX) % CONST::ANX;
			for (int k = any - ii; k <= any + ii; k++) {
				aiy = (k + CONST::ANY) % CONST::ANY;
				for (int l = anz - ii; l <= anz + ii; l++) {
					aiz = (l + CONST::ANZ) % CONST::ANZ;
					int lat_num = area_cell_num[aix][aiy][aiz];
					for (int m = 0; m < lat_num; m++) {
#ifdef TEST

						//break;
#endif
						o = area_cell_idx[aix][aiy][aiz][m];
						if (i == o || (membflg&&env->cell_state[o] == MEMB))continue; //not contains disa
						//we can eliminate this if
						
						diffx = cellx - ret_x(cellx, env->cell_x[o]);
						diffy = celly - ret_y(celly, env->cell_y[o]);
						diffz = cellz - env->cell_z[o];
						if (CONST::LJ_THRESH*(cell_rad + env->cell_radius[o])*t_sqrtD(diffx*diffx + diffy*diffy + diffz*diffz) >= 1) {
							env->cell_connected_index[i][env->cell_connected_num[i]] = o;
							env->cell_connected_num[i]++;
							assert(env->cell_connected_num[i] < CONST::N2);
						}
						
					}
				}
			}
		}
	}


}

void update_cell_dyn(ENV* env, ENV* new_env) {
	static double inv_dn[CONST::NMX*CONST::NMY], inv_dm[CONST::NMX*CONST::NMY], 
		nx[CONST::NMX*CONST::NMY], mx[CONST::NMX*CONST::NMY],
		ny[CONST::NMX*CONST::NMY], my[CONST::NMX*CONST::NMY],
		nz[CONST::NMX*CONST::NMY], mz[CONST::NMX*CONST::NMY],
		ipn[CONST::NMX*CONST::NMY], ipm[CONST::NMX*CONST::NMY];
	
	bend(env,new_env, nx, ny, nz, mx, my, mz,
		inv_dn, inv_dm, ipn, ipm);
	
	memb_memb_interact(env, new_env);
	der_der_interact(env, new_env);
	other_interact(env, new_env);
	pair_interact(env, new_env);
	fixed_memb_pos(env, new_env);
	state_renew(env, new_env);
	pair_disperse(env, new_env);
	pos_fix_periodic_cond(new_env);
	connect_lj(new_env);


}
