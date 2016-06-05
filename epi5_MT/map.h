#pragma once
#include "define.h"
#include "Field.h"
class CellManager;
void setup_map_lat(CellManager& cman,
	Field<Cell*,cont::NX + 1, cont::NY + 1, cont::NZ + 1>& cmap1,
	Field<uint_fast8_t,cont::NX + 1, cont::NY + 1, cont::NZ + 1>& cmap2 );