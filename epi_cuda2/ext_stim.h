#pragma once
#include "define.h"
class CellManager;
void calc_ext_stim(CellManager*cm, real* c_ext_stim, real* new_ext_stim, int* cmap1, float* cmap2, real zmax);