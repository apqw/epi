#pragma once
#include "define.h"
double calc_cell_P(const ENV& current_env, int update_cell);

void update_cell_internal(const ENV& current_env, ENV& updated_env);
void update_ca(const ENV* current_env,ENV* updated_env);
void update_map(ENV* update_env);
void set_air_stim(ENV* current_env);
void reset_ca_ave_in_cell(ENV* update_env);

void update_all(ENV** current_env,ENV** updated_env);
void cell_div(ENV* env,int cell_index);

void set_memv_indices(ENV* env);
void bend(ENV* env);
void bend_interac(ENV* env);
void memb_memb_interact(const ENV* env,ENV* new_env);
void der_der_interact(const ENV* env, ENV* new_env);
void other_interact(const ENV* env, ENV* new_env);
