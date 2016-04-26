#pragma once
#include "define.h"
double calc_cell_P(const ENV& current_env, int update_cell);

void update_cell_internal(const ENV& current_env, ENV& updated_env);
void update_ca(const ENV* current_env,ENV* updated_env);
void update_map(ENV* update_env);
void set_air_stim(ENV* current_env);
void reset_ca_ave_in_cell(ENV* update_env);

void update_all(ENV** current_env,ENV** updated_env);
void cell_div(ENV* env,ENV *update_env, int cell_idx, int given_divtimes, int dermis_index);

void set_memv_indices(ENV* env);
void bend(const ENV* env,ENV* new_env, double nx[], double ny[], double nz[], double mx[], double my[], double mz[],
	double inv_dn[], double inv_dm[], double ipn[], double ipm[]);
void bend_interac(ENV* env);
void memb_memb_interact(const ENV* env,ENV* new_env);
void der_der_interact(const ENV* env, ENV* new_env);
void other_interact(const ENV* env, ENV* new_env);
void pair_interact(const ENV* env, ENV* new_env);
void state_renew(ENV* env, ENV* update_env);
void pair_disperse(const ENV* env, ENV* new_env);
void pos_fix_periodic_cond(ENV* env);
void connect_lj(ENV* env);

void update_cell_dyn(ENV* env, ENV* new_env);
