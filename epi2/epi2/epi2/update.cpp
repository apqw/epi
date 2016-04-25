//#include <stdlib.h>
#include <stdio.h>
#include "update.h"
#include <math.h>
#include <omp.h>
#include "SFMT.h"
#define max(a,b) ((a)<=(b)?(b):(a))
#define min(a,b) ((a)<=(b)?(a):(b))

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
			if (env.cell_state[j] == ALIVE || env.cell_state[j] == FIX || env.cell_state[j] == MUSUME) {
				double ix = env.cell_x_index[j], iy = env.cell_y_index[j], iz = env.cell_z_index[j];
				double avg_Ca = calc_avg8(env.Ca2P_value, ix, iy, iz);
				double avg_B = calc_avg8(env.B_value, ix, iy, iz);
				double age = (env.cell_state[j] == FIX || env.cell_state[j] == MUSUME) ? env.cell_ageb[j] : env.cell_agek[j];
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
				if (env.cell_state[connected_cell_index] == ALIVE
					|| env.cell_state[connected_cell_index] == DEAD
					|| env.cell_state[connected_cell_index] == FIX
					|| env.cell_state[connected_cell_index] == MUSUME) {
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

void update_ca(const ENV *current_env,ENV *updated_env){
    int iz_bound=(int)((current_env->zzmax+CONST::FAC_MAP*CONST::R_max)*CONST::inv_dz);
    #pragma omp parallel for
    for(int j=0;j<CONST::NX;j++){
        __m256d mul1,mul2,sub1_1,sub1_2,sub2;
        alignas(32) double diffuse_coef[4];
        int prev_x = (j-1 +CONST::NX)%CONST::NX;
        int next_x = (j+1 + CONST::NX)%CONST::NX;
        for(int k=0;k<CONST::NY;k++){
            int prev_y = (k-1+CONST::NY)%CONST::NY;
            int next_y = (k+1+CONST::NY)%CONST::NY;
            for(int l=0;l<iz_bound;l++){
                int prev_z = l==0?1:l-1;
                int next_z = l==CONST::NZ?CONST::NZ-1:l+1;
                int map = current_env->cell_map[j][k][l];
                double diffu=0;
                if(map!=-1){
                    if(current_env->cell_state[map] == ALIVE||current_env->cell_state[map] == FIX||current_env->cell_state[map] == MUSUME){
                      diffu=current_env->cell_diff_c[map];
                    }
                }
               mul1 = _mm256_set_pd(
                            current_env->cell_map2[prev_x][k][l],
                            current_env->cell_map2[next_x][k][l],
                            current_env->cell_map2[j][prev_y][l],
                            current_env->cell_map2[j][next_y][l]);
                mul2 = _mm256_set_pd(
                            current_env->cell_map2[j][k][prev_z],
                            current_env->cell_map2[j][k][next_z],
                            0,
                            0
                            );
                sub1_1 = _mm256_set_pd(
                            current_env->Ca2P_value[prev_x][k][l],
                            current_env->Ca2P_value[next_x][k][l],
                            current_env->Ca2P_value[j][prev_y][l],
                            current_env->Ca2P_value[j][next_y][l]
                            );
                sub1_2 = _mm256_set_pd(
                            current_env->Ca2P_value[j][k][prev_z],
                            current_env->Ca2P_value[j][k][next_z],
                            0,
                            0
                            );
                sub2 = _mm256_set1_pd(current_env->Ca2P_value[j][k][l]);

                _mm256_store_pd(diffuse_coef,_mm256_add_pd(_mm256_mul_pd(mul1,_mm256_sub_pd(sub1_1,sub2)),_mm256_mul_pd(mul2,_mm256_sub_pd(sub1_2,sub2))));
                updated_env->Ca2P_value[j][k][l] = current_env->Ca2P_value[j][k][l]
                        +CONST::dt_ca*(CONST::Ca2P::dA*(
                                           diffuse_coef[0]+diffuse_coef[1]+diffuse_coef[2]+diffuse_coef[3])
                        *CONST::inv_dxSq+CONST::Ca2P::Kac*diffu-CONST::Ca2P::Kaa*current_env->Ca2P_value[j][k][l]+current_env->air_stim[j][k][l]);
                /*
                updated_env->Ca2P_value[j][k][l] = current_env->Ca2P_value[j][k][l]
                        +CONST::dt_ca*(CONST::Ca2P::dA*(
                                           current_env->cell_map2[prev_x][k][l]*(current_env->Ca2P_value[prev_x][k][l]-current_env->Ca2P_value[j][k][l])
                                           +current_env->cell_map2[next_x][k][l]*(current_env->Ca2P_value[next_x][k][l]-current_env->Ca2P_value[j][k][l])
                                           +current_env->cell_map2[j][prev_y][l]*(current_env->Ca2P_value[j][prev_y][l]-current_env->Ca2P_value[j][k][l])
                                           +current_env->cell_map2[j][next_y][l]*(current_env->Ca2P_value[j][next_y][l]-current_env->Ca2P_value[j][k][l])
                                           +current_env->cell_map2[j][k][prev_z]*(current_env->Ca2P_value[j][k][prev_z]-current_env->Ca2P_value[j][k][l])
                                           +current_env->cell_map2[j][k][next_z]*(current_env->Ca2P_value[j][k][next_z]-current_env->Ca2P_value[j][k][l])
                                           )*CONST::inv_dxSq+CONST::Ca2P::Kac*diffu-CONST::Ca2P::Kaa*current_env->Ca2P_value[j][k][l]);
                                           */


            }
        }
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
                            update_env->cell_map[ipx][ipy][iz+j]=i;
                        }
                    }
                }
            }
            if(update_env->cell_state[i] ==FIX||update_env->cell_state[i] ==MUSUME||update_env->cell_state[i] ==ALIVE){
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
                                update_env->cell_map2[ipx][ipy][iz+j]=1;
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
    update_map(*current_env);
    reset_ca_ave_in_cell(*updated_env); //do not use current value
    set_air_stim(*current_env);
    ENV* old_env,*new_env,*tmp_env_ptr;
    old_env=*current_env;new_env=*updated_env;
    for (int i = 0; i < CONST::Ca2P::ca_N; i++) {
        printf("%d\n", i);
        update_cell_internal(*old_env,*new_env);
        update_ca(old_env,new_env);

        tmp_env_ptr=new_env;
        new_env=old_env;
        old_env=tmp_env_ptr; //old is latest now
    }
    for(int i=old_env->NMEMB+old_env->NDER;i<old_env->current_cell_num;i++){
        if(old_env->cell_state[i] == ALIVE || old_env->cell_state[i] == FIX || old_env->cell_state[i] == MUSUME){
            old_env->cell_ca_ave[i]/=CONST::Ca2P::ca_N;
        }else{
            old_env->cell_ca_ave[i]=old_env->cell_c[i]=0.0;
        }
    }

}

void reset_ca_ave_in_cell(ENV *update_env){
    for(int i=update_env->NMEMB+update_env->NDER;i<update_env->current_cell_num;i++){
        if(update_env->cell_state[i]==ALIVE||
                update_env->cell_state[i]==FIX||
                update_env->cell_state[i]==MUSUME){
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

void bend(ENV* env,double nx[],double ny[],double nz[], double mx[], double my[], double mz[],
	double inv_dn[],double inv_dm[],double ipn[],double ipm[]) {
	for (int j = 0; j < env->NMEMB; j++) {
		
		double diffx, diffy, diffz;
		diffx = env->cell_x[env->memb_idx_jr[j]] - env->cell_x[j];
		diffy = env->cell_y[env->memb_idx_jr[j]] - env->cell_y[j];
		diffz = env->cell_z[env->memb_idx_jr[j]] - env->cell_z[j];

		//replace this real sqrt func in production env
		inv_dn[j] = t_sqrtD(diffx*diffx + diffy*diffy + diffz*diffz);
		if (inv_dn[j] > 1.0e+5) {
			printf("error: distance dn too short\n");
			exit(1);
		}
		nx[j] = perio_diff(env->cell_x[env->memb_idx_jr[j]], env->cell_x[j], CONST::Lx)*inv_dn[j];
		ny[j] = perio_diff(env->cell_y[env->memb_idx_jr[j]], env->cell_y[j], CONST::Ly)*inv_dn[j];
		nz[j] = diffz*inv_dn[j];

		diffx = env->cell_x[env->memb_idx_ju[j]] - env->cell_x[j];
		diffy = env->cell_y[env->memb_idx_ju[j]] - env->cell_y[j];
		diffz = env->cell_z[env->memb_idx_ju[j]] - env->cell_z[j];
		inv_dm[j] = t_sqrtD(diffx*diffx + diffy*diffy + diffz*diffz);
		if (inv_dm[j] > 1.0e+5) {
			printf("error: distance dn too short\n");
			exit(1);
		}
		mx[j] = perio_diff(env->cell_x[env->memb_idx_ju[j]], env->cell_x[j], CONST::Lx)*inv_dm[j];
		my[j] = perio_diff(env->cell_y[env->memb_idx_ju[j]], env->cell_y[j], CONST::Ly)*inv_dm[j];
		mz[j] = diffz*inv_dm[j];

	}
	for (int j = 0; j < env->NMEMB; j++) {
		ipn[j] = nx[env->memb_idx_jr[j]] * nx[j] + ny[env->memb_idx_jr[j]] * ny[j] + nz[env->memb_idx_jr[j]] * nz[j];
		ipm[j] = mx[env->memb_idx_ju[j]] * mx[j] + my[env->memb_idx_ju[j]] * my[j] + mz[env->memb_idx_ju[j]] * mz[j];
	}

	for (int j = 0; j < env->NMEMB; j++) {
		env->cell_x[j] += CONST::dt_cell*bend_force_sqr(nx, mx, ipn, ipm, inv_dn, inv_dm, j, env->memb_idx_jr[j], env->memb_idx_jl[j], env->memb_idx_jll[j], env->memb_idx_ju[j], env->memb_idx_jb[j], env->memb_idx_jbb[j]);
		env->cell_y[j] += CONST::dt_cell*bend_force_sqr(ny, my, ipn, ipm, inv_dn, inv_dm, j, env->memb_idx_jr[j], env->memb_idx_jl[j], env->memb_idx_jll[j], env->memb_idx_ju[j], env->memb_idx_jb[j], env->memb_idx_jbb[j]);
		env->cell_z[j] += CONST::dt_cell*bend_force_sqr(nz, mz, ipn, ipm, inv_dn, inv_dm, j, env->memb_idx_jr[j], env->memb_idx_jl[j], env->memb_idx_jll[j], env->memb_idx_ju[j], env->memb_idx_jb[j], env->memb_idx_jbb[j]);
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

inline double ret__(double x1, double x2,double axis_len)
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

void memb_memb_interact(const ENV* env, ENV* new_env) {
	for (int j = 0; j < env->NMEMB; j++) {
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

			double diffx = env->cell_x[j] - env->cell_x[connected_index];
			double diffy = env->cell_y[j] - env->cell_y[connected_index];
			double diffz = env->cell_z[j] - env->cell_z[connected_index];

			double distlj_inv = t_sqrtD(diffx*diffx + diffy*diffy+diffz*diffz);
			double rad_sum = (env->cell_radius[j] + env->cell_radius[connected_index]);
			double cr_dist = rad_sum*CONST::P_MEMB;
			double cr_dist_inv = 1 / cr_dist;
			double lambda_dist = (1.0 + CONST::P_MEMB)*rad_sum;
			double cr_lj = cr_dist*distlj_inv;
			double dum_lj6,dum_ljmain;
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
			double dumxx = CONST::dt_cell*dum_ljmain*(env->cell_x[j] - ret__(env->cell_x[j], env->cell_x[connected_index],CONST::Lx));
			double dumyy = CONST::dt_cell*dum_ljmain*(env->cell_y[j] - ret__(env->cell_y[j], env->cell_y[connected_index], CONST::Ly));
			double dumzz = CONST::dt_cell*dum_ljmain*(env->cell_z[j] - env->cell_z[connected_index]);
			new_env->cell_x[j] += dumxx;
			new_env->cell_y[j] += dumyy;
			new_env->cell_z[j] += dumzz;
			new_env->cell_x[connected_index] -= dumxx;
			new_env->cell_y[connected_index] -= dumyy;
			new_env->cell_z[connected_index] -= dumzz;
		}
	}
}
//celldyn
void cell_div(ENV *env,int target_idx){
	genrand_res53();
    if(env->cell_div_pair_index[target_idx]!=-1)return;
    //double par=env->cell_state[target_idx]==MUSUME?CONST::ag
}
