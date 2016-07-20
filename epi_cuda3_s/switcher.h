#pragma once
#include "CData.h"
#include <type_traits>
class switcher{
	CValue<int> _swt;
public:
	__host__ void switch_value();
	__host__ __device__ int current()const;
	__host__ __device__ int next()const;
	switcher();
};

template<class T>
class SwtAccessor{
	const switcher& _swt;
	T& ref;
public:
	__device__ __host__ auto current()->decltype(ref[_swt.current()]){
		return ref[_swt.current()];
	}
	__device__ __host__ auto current() const->typename std::add_const<decltype(ref[_swt.current()])>::type {
		return ref[_swt.current()];
	}

	__device__ __host__ auto next()->decltype(ref[_swt.next()]){
		return ref[_swt.next()];
	}

	__device__ __host__ auto next() const->typename std::add_const<decltype(ref[_swt.next()])>::type {
		return ref[_swt.current()];
	}

	template<class F>
	__device__ __host__ void foreach(F&& lmbd){
		lmbd(current());
		lmbd(next());
	}
	template<typename U>
	__device__ __host__ void set(size_t index, U&& item){
		current()[index] = item;
		next()[index] = item;
	}
	SwtAccessor(const switcher& swt,T& data_ref) :_swt(swt), ref(data_ref){}

};