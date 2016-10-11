
#ifndef MAP_H_
#define MAP_H_
#include <cstdint>
#include "../../misc/DynArr.h"
class CellManager;
class Cell;
void setup_map_lat(CellManager & cman, Dyn3DArr<const Cell*>& cmap1, Dyn3DArr<uint_fast8_t>& cmap2);

#endif /* MAP_H_ */
