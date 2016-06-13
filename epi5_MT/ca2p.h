#pragma once
#include "swapdata.h"
#include "Field.h"
#include "define.h"
#include "cellmanager.h"

void calc_ca2p(CellManager& cman, SwapData<FArr3D<double>>& ATP, const FArr3D<double>& ext_stim_first,const FArr3D<const Cell*>& cmap1, const FArr3D<cmask_ty>& cmap2, double zzmax);