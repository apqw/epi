#pragma once
#include <memory>
#include "atomics.h"
#include "define.h"
class Cell;
using CellPtr = std::shared_ptr<Cell>;
class CellManager :public Lockfree_push_stack<CellPtr, cont::MAX_CELL_NUM> {
	Lockfree_push_stack<size_t, cont::MAX_CELL_NUM> remove_queue;
public:
	size_t register_cell(const CellPtr& c);
	void add_remove_queue(size_t idx);
	void remove_exec();
};