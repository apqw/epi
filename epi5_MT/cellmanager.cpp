#include "cellmanager.h"
#include "cell.h"

size_t CellManager::register_cell(const CellPtr & c)
{
	return push_back_with_index(c);
}

void CellManager::add_remove_queue(size_t idx)
{
	remove_queue.push_back(idx);
}

void CellManager::remove_exec()
{
	for (size_t i = 0; i < remove_queue.size(); ++i) {
		_data[i] = _data[--_next];
		_data[i]->migrate(i);
	}
	remove_queue.clear();
}
