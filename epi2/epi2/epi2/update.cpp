#include "update.h"

inline double calc_avg8(const Arr3D<double, CONST::NX, CONST::NY, CONST::NZ>& arr, int ix, int iy, int iz) {
	return 0.125*(arr[ix][iy][iz] + arr[ix][iy][iz + 1] +
		+arr[ix][iy + 1][iz] + arr[ix][iy + 1][iz + 1]
		+ arr[ix + 1][iy][iz] + arr[ix + 1][iy][iz + 1]
		+ arr[ix + 1][iy + 1][iz] + arr[ix + 1][iy + 1][iz + 1]);
}

inline double fu(double u, double v, double p, double B) {
	//kgÇ∆KgÇ™Ç†ÇÈÅ@âΩåÃ
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

	for (int i = 0; i < CONST::Ca2P::ca_N; i++) {
		double diffu[CONST::MAX_CELL_NUM] = { 0 };
		printf("%d\n", i);
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
			double w_diff = 0;
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
			updated_env.cell_P[j] = env.cell_P[j] + CONST::dt_ca*(p_diff + CONST::Ca2P::dP*IAG_val*p_diff_unmult);
			updated_env.cell_h[j] = env.cell_h[j] + CONST::dt_ca*v_diff;
			

		}


	}
}