#ifndef CELL_H
#define CELL_H
#include "define.h"
#include <memory>
#include <vector>
#include <tbb/concurrent_vector.h>

class Cell
{
public:
    CELL_STATE state;
    double age;
    double pos[3];
    double nat_spr_len;
    Cell();
};
using CellPtr = std::shared_ptr<Cell>;
class CellMan{
tbb::concurrent_vector<CellPtr> cell;
public:
void add(CellPtr& c);
void remove(int index);

};

#endif // CELL_H
