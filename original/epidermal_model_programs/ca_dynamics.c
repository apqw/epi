#include <stdlib.h>
#include <math.h>
#include "param.h"
#include "funcs.h"

void ca_dynamics(double *u_ave, double B[NX + 1][NY + 1][NZ + 1],
	double *u, double *v, double *p, double a[NX + 1][NY + 1][NZ + 1],
	double w[NN][N2], int *ixx, int *iyy, int *izz, int lj[NN][N2], int indx[],
	int *state, double *ageb, double *agek, int map[NX + 1][NY + 1][NZ + 1],
	int map2[NX + 1][NY + 1][NZ + 1], double zzmax,
	double *r, double *xx, double *yy, double *zz, int n)
{
	int j, k, l, m, ii;
	int ix, iy, iz;
	int prev_x, prev_y, prev_z, next_x, next_y, next_z;
	double tmp_a, tmp_B, diff, tmp;
	double old_u[NN], old_v[NN], old_p[NN], old_fat[NN], old_Vc[NN], old_w[NN][N2], old_a[NX + 1][NY + 1][NZ + 1];
	int iz_bound = (int)((zzmax + FAC_MAP*R_max) / dz);

	double diffu[NN];
	double dbg;
	double air_stim[NX + 1][NY + 1][NZ + 1];

	for (j = NMEMB + NDER; j < n; j++)
		if (state[j] == ALIVE || state[j] == FIX || state[j] == MUSUME)
			u_ave[j] = 0.0;

	// air exposure
	setup_air_stim(map, lj, indx, state, air_stim);

	for (ii = 0; ii < ca_N; ii++){

		for (j = NMEMB + NDER; j < n; j++) {
			old_u[j] = u[j];
			old_v[j] = v[j];
			old_p[j] = p[j];
			diffu[j] = 0.0;

			for (k = 0; k < N2; k++) {
				old_w[j][k] = w[j][k];
			}
		}
		for (j = 0; j <= NX; j++) {
			for (k = 0; k <= NY; k++) {
				for (l = 0; l <= iz_bound; l++) {
					old_a[j][k][l] = a[j][k][l];
				}
			}
		}

		for (j = NMEMB + NDER; j < n; j++) {
			if (state[j] == DEAD) {
				p[j] += -DT * J(old_p[j]);//me_dead_val
				for (k = 0; k < indx[j]; k++) {
					if ((l = lj[j][k]) > -1 && (state[l] == ALIVE))
						p[j] += DT * dp * (old_p[l] - old_p[j]);//react_alive_val
				}
				//まずold_p[l]は普通に足す、足した回数を数えてその回数*old_p[j]を引く、そのあとdpを賭ける
				//このループを抜けた後にDTをかけるならかける、一番最後に賭けても良い気がする
			}
		}

		for (j = NMEMB + NDER; j < n; j++) {
			if (state[j] == ALIVE || state[j] == FIX || state[j] == MUSUME) {

				ix = ixx[j];
				iy = iyy[j];
				iz = izz[j];

				if (iz >= iz_bound) {
					printf("error: illegal iz value in function ca_dynamics.c\n");
					exit(1);
				}

				tmp_a = (old_a[ix][iy][iz] + old_a[ix + 1][iy][iz] + old_a[ix][iy + 1][iz]
					+ old_a[ix][iy][iz + 1] + old_a[ix + 1][iy + 1][iz] + old_a[ix + 1][iy][iz + 1]
					+ old_a[ix][iy + 1][iz + 1] + old_a[ix + 1][iy + 1][iz + 1]) / 8.0;
				tmp_B = (B[ix][iy][iz] + B[ix + 1][iy][iz] + B[ix][iy + 1][iz]
					+ B[ix][iy][iz + 1] + B[ix + 1][iy + 1][iz] + B[ix + 1][iy][iz + 1]
					+ B[ix][iy + 1][iz + 1] + B[ix + 1][iy + 1][iz + 1]) / 8.0;

				diffu[j] = fu(old_u[j], old_v[j], old_p[j], tmp_a, tmp_B, j, ii);

				v[j] += DT * fv(old_u[j], old_v[j], old_p[j], tmp_a, ((state[j] == FIX || state[j] == MUSUME) ? ageb[j] : agek[j]), state[j]);
				p[j] += DT * fp(old_u[j], old_v[j], old_p[j], tmp_a, ((state[j] == FIX || state[j] == MUSUME) ? ageb[j] : agek[j]), state[j]);
				//fp 後ろ三つしか使ってない
				//me_val
				for (k = 0; k < indx[j]; k++) {
					l = lj[j][k];
					if ((state[l] == ALIVE || state[l] == DEAD || state[l] == FIX || state[l] == MUSUME)) {

						diffu[j] += du * IAG(agek[j], state[j]) * old_w[j][k] * (old_u[l] - old_u[j]);
						p[j] += DT * dp * IAG(agek[j], state[j]) * old_w[j][k] * (old_p[l] - old_p[j]); //react_val

						if (state[l] == ALIVE) {
							diff = fabs(old_u[l] - old_u[j]);
							w[j][k] += DT * fw(diff, old_w[j][k]);
						}
					}
				}
				u[j] += DT * diffu[j];
				u_ave[j] += u[j];
			}

		}
		/*ATP計算*/
		for (j = 0; j < NX; j++){
			for (k = 0; k < NY; k++) {
				//periodic
				prev_x = (j - 1 + NX) % NX;
				next_x = (j + 1 + NX) % NX;
				prev_y = (k - 1 + NY) % NY;
				next_y = (k + 1 + NY) % NY;
				for (l = 0; l < iz_bound; l++) {
					//neumann?
					if (l == 0) {
						prev_z = 1;
					}
					else {
						prev_z = l - 1;
					}
					if (l == NZ) {
						next_z = NZ - 1;
					}
					else {
						next_z = l + 1;
					}
					m = map[j][k][l];
					if (m == -1) tmp = 0.0;
					else {
						if (state[m] == ALIVE || state[m] == FIX || state[m] == MUSUME)
							tmp = diffu[m];
						else tmp = 0.0;
					}
					//実装おかしい？
					//Sum[i=1,N]G~~の部分は？
					a[j][k][l] = old_a[j][k][l]
						+ DT * (Da * (map2[prev_x][k][l] * (old_a[prev_x][k][l] - old_a[j][k][l])
						+ map2[next_x][k][l] * (old_a[next_x][k][l] - old_a[j][k][l])
						+ map2[j][prev_y][l] * (old_a[j][prev_y][l] - old_a[j][k][l])
						+ map2[j][next_y][l] * (old_a[j][next_y][l] - old_a[j][k][l])
						+ map2[j][k][prev_z] * (old_a[j][k][prev_z] - old_a[j][k][l])
						+ map2[j][k][next_z] * (old_a[j][k][next_z] - old_a[j][k][l])) / (dx*dx)
						+ fa(tmp, old_a[j][k][l]));
					// air exposure
					a[j][k][l] += DT * air_stim[j][k][l];
				}
			}
		}
		/* B.C. */
		for (l = 0; l <= iz_bound; l++) {
			for (j = 0; j < NX; j++) a[j][NY][l] = a[j][0][l];
			for (k = 0; k <= NY; k++) a[NX][k][l] = a[0][k][l];
		}

	}

	for (j = NMEMB + NDER; j < n; j++) {
		if (state[j] == ALIVE || state[j] == FIX || state[j] == MUSUME)
			u_ave[j] /= ca_N;
		else
			u_ave[j] = u[j] = 0.0;
	}

}

double J(double p)
{
	return para_kpp * p;
}

/* reaction terms for ATP */
double fa(double diffu, double a)
{
	return I(diffu) - g(a);
}
double I(double diffu)
{
	double dum = positive(diffu);
	return STIM11 * dum;
}


double g(double a)
{
	return para_kaa * a;
}

/* reaction terms for IP3 */
double fp(double u, double v, double p, double a, double age, int state)
{
	double h(double, double, int);
	return h(a, age, state) - J(p);
}
double h(double a, double age, int state)
{
	double Kpa(double, int);
	return Kpa(age, state) * a / (para_H0 + a);
}
double Kpa(double age, int state)
{
	if (state == ALIVE){
		return para_Kgra + 0.5*(para_Kpri - para_Kgra) * (1.0 + tanh((THRESH_SP - age) / delta_K));
	}
	else if (state == FIX || state == MUSUME){
		return para_Kpri;
	}
	else {
		return 0;
	}
}

/* reaction terms for Ca */
double fu(double u, double v, double p, double a, double B, int debug_j, int ii)
{
	//  static double result_debug=0.0;

	return para_kf * (mu0 + mu1 * p / (p + para_kmu)) * (para_b + para_bb * u / (para_k1 + u)) * v
		- gamma * u / (para_kg + u) + beta + para_kbc * B*B * para_Cout / (para_Hb + B*B);

}

/* reaction terms for h */
double fv(double u, double v, double p, double a, double age, int state)
{
	double th(double, int);
	return (para_k2 * para_k2 / (para_k2 * para_k2 + u * u) - v) / th(age, state);
}
double th(double age, int state)
{
	if (state == ALIVE){
		return para_thgra + ((para_thpri - para_thgra) / 2.) * (1.0 + tanh((THRESH_SP - age) / delta_th));
	}
	else if (state == FIX || state == MUSUME){
		return para_thpri;
	}
	else {
		printf("error in function fv\n");
		exit(1);
	}
}

/*GJの閉じ具合調節*/
double fw(double diff, double w)
{
	return (1. - w) + (-1. + tanh((para_wd - diff) / 0.1)) / 2.; //epsw0 == 0.1 //done
	//-1. <-????
}

/* inter-cellular exchange */
double IAG(double age, int state)
{
	if (state == FIX || state == MUSUME){
		return para_iage_kitei;
	}
	else {
		return 0.5*(1.0 + tanh((age - THRESH_SP) / delta_I));
	}
}

void setup_air_stim(int map[NX + 1][NY + 1][NZ + 1], int lj[][N2], int indx[], int state[], double air_stim[NX + 1][NY + 1][NZ + 1])
{
	int j, k, l, jj, kk, ll;

	for (j = 0; j <= NX; j++)
		for (k = 0; k <= NY; k++)
			for (l = 0; l <= NZ; l++) {
				air_stim[j][k][l] = 0.0;
				jj = map[j][k][l];
				if (jj != -1) {
					for (kk = 0; kk < indx[jj]; kk++) {
						if (state[lj[jj][kk]] == AIR) {
							air_stim[j][k][l] = AIR_STIM;
							break;
						}
					}
				}
			}
}

