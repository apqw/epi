#pragma once
#include "global_cuda.h"
class remove_queue
{
	
public:
	unsigned int* _atomic_queue_head;
	int* queue;
	unsigned int* queue_max_size;
	remove_queue(unsigned int queue_size);
	~remove_queue();
	__device__ void push_back(int index);
	__device__ unsigned int queue_head()const;
	__device__ void reset();
};

