#pragma once
#include "define.h"
double calc_cell_P(const ENV& current_env, int update_cell);

void update_cell_internal(const ENV& current_env, ENV& updated_env);