//#include <stdlib.h>
#include <stdio.h>
#include "update.h"
#include <math.h>
#include <omp.h>
#include "SFMT.h"
#define max(a,b) ((a)<=(b)?(b):(a))
#define min(a,b) ((a)<=(b)?(a):(b))
inline double calc_avg8(const double arr[CONST::NX+1][CONST::NY+1][CONST::NZ+1], int ix, int iy, int iz) {
	return 0.125*(arr[ix][iy][iz] + arr[ix][iy][iz + 1] +
		+arr[ix][iy + 1][iz] + arr[ix][iy + 1][iz + 1]
		+ arr[ix + 1][iy][iz] + arr[ix + 1][iy][iz + 1]
		+ arr[ix + 1][iy + 1][iz] + arr[ix + 1][iy + 1][iz + 1]);
}

inline double fu(double u, double v, double p, double B) {
	//kg‚ÆKg‚ª‚ ‚é@‰½ŒÌ
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

//celldyn
void cell_div(ENV *env,int target_idx){
    if(env->cell_div_pair_index[target_idx]!=-1)return;
    //double par=env->cell_state[target_idx]==MUSUME?CONST::ag
}
