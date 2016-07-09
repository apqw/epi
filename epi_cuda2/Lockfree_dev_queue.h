#pragma once
#include "define.h"
#include <device_functions.h>
//standard layout
template<typename T,size_t N>
class LockfreeDeviceQueue{
public:
	unsigned int head;

	T data[N];
	LockfreeDeviceQueue() :head(0){}

	__device__ void push_back(T item);
	__device__ unsigned int size()const{
		return head;
	}
	__device__ void reset(){
		head = 0;
	}
};

