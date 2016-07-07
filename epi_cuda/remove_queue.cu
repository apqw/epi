#include "remove_queue.cuh"

remove_queue::remove_queue(unsigned int queue_size)
{
	cudaMalloc((void**)&queue, sizeof(int)*queue_size);
	cudaMalloc((void**)&_atomic_queue_head, sizeof(unsigned int));
	cudaMemset(_atomic_queue_head, 0, sizeof(unsigned int));
	cudaMalloc((void**)&queue_max_size, sizeof(unsigned int));
	cudaMemcpy(queue_max_size, &queue_size, sizeof(unsigned int),cudaMemcpyHostToDevice);
}

remove_queue::~remove_queue(){
	cudaFree(queue);
	cudaFree(_atomic_queue_head);
	cudaFree(queue_max_size);
}

__device__ void remove_queue::push_back(int index){
	queue[atomicAdd(_atomic_queue_head, 1)] = index;
}

__device__ unsigned int remove_queue::queue_head()const{
	return *_atomic_queue_head;
}

__device__ void remove_queue::reset(){
	*_atomic_queue_head = 0;
}