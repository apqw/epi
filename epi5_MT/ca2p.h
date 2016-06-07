#pragma once
#include "swapdata.h"
#include "Field.h"
#include "define.h"
#include "cellmanager.h"

void calc_ca2p(CellManager& cman, SwapData<FArr3D<double>>& ca2p, const  SwapData<FArr3D<double>>& ext_stim, FArr3D<Cell*>& cmap1, FArr3D<uint_fast8_t>& cmap2, double zzmax);