#pragma once
#include "define.h"
#include "Field.h"
class CellManager;
void setup_map_lat(CellManager& cman,
	FArr3D<const Cell*>& cmap1,
	FArr3D<cmask_ty>& cmap2 );