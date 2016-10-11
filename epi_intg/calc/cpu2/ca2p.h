#pragma once
#include "../../misc/DynArr.h"
#include "../../misc/swapdata.h"
#include "../../define.h"
class CellManager;
class Cell;

void calc_ca2p(CellManager& cman, SwapData<Dyn3DArr<real>>& ATP, const Dyn3DArr<real>& ext_stim_first, const Dyn3DArr<const Cell*>& cmap1, const Dyn3DArr<uint_fast8_t>& cmap2, real zzmax);
