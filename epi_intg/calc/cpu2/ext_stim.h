#pragma once
#ifndef EXT_STIM_H_
#define EXT_STIM_H_
#include "../../misc/swapdata.h"
#include "../../misc/DynArr.h"
#include "../../define.h"
class Cell;
void calc_ext_stim(SwapData<Dyn3DArr<real>>& ext_stim, const Dyn3DArr<const Cell*>& cmap1, const Dyn3DArr<uint_fast8_t>& cmap2, real zzmax);
#endif