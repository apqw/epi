/*
 * utilc.cpp
 *
 *  Created on: 2016/07/09
 *      Author: yasu7890v
 */




#include "utils.h"

std::string int_to_string(int number)
{
  std::stringstream ss;
  ss << number;
  return ss.str();
}

#if __CUDACC_VER_MAJOR__ < 8 && !defined(USE_FLOAT)
__device__ double atomicAdd(double* address, double val)
{
	unsigned long long int* address_as_ull =
		(unsigned long long int*)address;
	unsigned long long int old = *address_as_ull, assumed;
	do {
		assumed = old;
		old = atomicCAS(address_as_ull, assumed,
			__double_as_longlong(val +
			__longlong_as_double(assumed)));
	} while (assumed != old);
	return __longlong_as_double(old);
}
#endif