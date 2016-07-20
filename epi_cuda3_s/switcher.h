#pragma once
#include "CData.h"

class switcher{
	CValue<int> _swt;
public:
	__host__ void switch_value();
	__host__ __device__ int current()const;
	__host__ __device__ int next()const;
	switcher();
};