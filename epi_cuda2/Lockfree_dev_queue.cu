#include "Lockfree_dev_queue.h"
template class LockfreeDeviceQueue<CellIndex, MAX_CELL_NUM>;
template<typename T,size_t N>
__device__ void LockfreeDeviceQueue<T,N>::push_back(T item){
	data[atomicAdd(&head, 1)] = item;
}