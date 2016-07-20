#pragma once
#include "define.h"
#include <thrust/device_malloc.h>
#include <thrust/device_free.h>
#include <thrust/fill.h>
#include <vector>
#include <memory>

template<typename T, size_t N = MAX_CELL_NUM>
struct CArr{
	devPtr<T> ptr;
	const bool no_free;
	__device__ __host__ operator devPtr<T>(){
		return ptr;
	}
	__device__ __host__ devRef<T> operator[](size_t index){
		return ptr[index];
	}
	CArr() :ptr(thrust::device_malloc<T>(N)),no_free(false){
#ifdef DBG
		printf("CArr default ctor called.\n");
#endif
	}

	CArr(devPtr<T>& _ptr,bool _no_free=false) :ptr(_ptr),no_free(_no_free){
#ifdef DBG
		printf("CArr ctor with ptr called.\n");
#endif
	}
	CArr(devPtr<T>&& _ptr, bool _no_free = false) :ptr(_ptr), no_free(_no_free){
#ifdef DBG
		printf("CArr ctor with rvalue ptr called.\n");
#endif
	}

	~CArr(){
		if (!no_free){
#ifdef DBG
			printf("CArr dtor called and freed.\n");
#endif
			thrust::device_free(ptr);
		}
		else{
#ifdef DBG
			printf("CArr dtor called and not freed.\n");
#endif
		}
	}

	__host__ void fill(const T& value){
		thrust::fill(ptr, ptr + N,value);
	}
	
};

template<typename T>
struct CValue :public CArr<T, 1>{
	__device__ __host__ operator T()const{
		return (T)(this->ptr[0]);
	}
	__device__ __host__ CValue& operator=(T&& item){
		this->ptr[0] = item;
		return *this;
	}
	CValue() :CArr<T,1>(){
#ifdef DBG
		printf("CValue default ctor called.\n");
#endif

	}
	CValue(const T& init){
#ifdef DBG
		printf("CValue ctor with initial value called.\n");
#endif
		this->ptr[0] = init;
	}
	
};

template<typename T, size_t X, size_t Y, size_t Z>
struct CArr3D :public CArr<T, X*Y*Z>{
	__device__ __host__ static size_t _midx(size_t _x, size_t _y, size_t _z){
		return _x*Y*Z + _y*Z + _z;
	}
	__device__ __host__ devRef<T> operator()(size_t x, size_t y, size_t z){
		return this->operator[](_midx(x, y, z));
	}
	CArr3D() :CArr<T,X*Y*Z>(){}
	CArr3D(devPtr<T>& _ptr, bool _no_free = false) :CArr<T,X*Y*Z>(_ptr, _no_free){}
	CArr3D(devPtr<T>&& _ptr, bool _no_free = false) :CArr<T,X*Y*Z>(_ptr, _no_free){}
};

template<typename T>
using Field3D = CArr3D<T, NX + 1, NY + 1, NZ + 1>;

template<typename T,size_t N,size_t Ch,typename Arr=CArr<T,N>>
struct CArrMulti{
private:
	std::vector<std::unique_ptr<Arr>> c;
	devPtr<T> ptr;
public:

	CArrMulti() :ptr(thrust::device_malloc<T>(N*Ch)){
		//c.reserve(Ch);
#ifdef DBG
		printf("CArrMulti ctor called.\n");
#endif
		for (size_t i = 0; i < Ch; i++){
			c.push_back(
					std::unique_ptr<Arr>(new Arr(ptr + i*N, true))
					);
		}
	}
	__host__ Arr& operator[](size_t ch){
#ifdef DBG
		assert(ch < Ch);
#endif
		return *c[ch];
	}

	~CArrMulti(){
#ifdef DBG
		printf("CArrMulti dtor called and freed.\n");
#endif
		thrust::device_free(ptr);
	}



};

template<typename T, size_t X, size_t Y, size_t Z,size_t Ch>
struct CArr3DMulti :public CArrMulti<T, X*Y*Z,Ch,CArr3D<T,X,Y,Z>>{
	
};
