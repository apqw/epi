#pragma once
#include "swapdata.h"
#include "Field.h"
#include "define.h"
void calc_ext_stim(SwapData<FArr3D<double>>& ext_stim,const FArr3D<const Cell*>& cmap1, const FArr3D<cmask_ty>& cmap2, double zzmax);