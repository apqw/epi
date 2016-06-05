#pragma once
#include "swapdata.h"
#include "Field.h"
#include "define.h"
void calc_ext_stim(SwapData<Field<double, cont::NX + 1, cont::NY + 1, cont::NZ + 1>>& ext_stim, Field<Cell*, cont::NX + 1, cont::NY + 1, cont::NZ + 1>& cmap1,
	Field<uint_fast8_t, cont::NX + 1, cont::NY + 1, cont::NZ + 1>& cmap2, double zzmax);