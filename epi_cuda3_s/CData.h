#pragma once
#include "define.h"
#include <thrust/device_malloc.h>
#include <thrust/device_free.h>
#include <thrust/fill.h>
#include <vector>
#include <memory>
#include <iostream>
#include "utils.h"
#include <array>
#include <typeinfo>
#include "fsys.h"
template<typename T, size_t N = MAX_CELL_NUM>
struct CArr{
protected:
	const bool no_free;
	bool initialized=false;
	devPtr<T> ptr;
public:
	typedef typename devPtr<T> ptr_type;
	typedef typename T value_type;

	__host__ bool has_initialized()const{
		return initialized;
	}

	__device__ __host__ operator devPtr<T>(){
		return ptr;
	}

	__device__ __host__ devRef<T> operator[](size_t index){
		return ptr[index];
	}
	__device__ __host__ const devRef<const T> operator[](size_t index)const{
		return ptr[index];
	}
	CArr() :ptr(thrust::device_malloc<T>(N)),no_free(false){
		dbgprint("CArr default ctor called.\n");
	}

	CArr(devPtr<T>& _ptr,bool _no_free=false) :ptr(_ptr),no_free(_no_free){
		dbgprint("CArr ctor with ptr called.\n");
	}
	CArr(devPtr<T>&& _ptr, bool _no_free = false) :ptr(_ptr), no_free(_no_free){
		dbgprint("CArr ctor with rvalue ptr called.\n");
	}

	~CArr(){
		if (!no_free){
			dbgprint("CArr dtor called and freed.\n");
			thrust::device_free(ptr);
		}
		else{
			dbgprint("CArr dtor called and not freed.\n");
		}
	}

	__host__ void fill(const T& value){
		thrust::fill(ptr, ptr + N,value);
	}
	template<typename F>
	__host__ void initialize(F&& lmbd){
		if (initialized){
			std::cout << "Already initialized. ADDR:" << (unsigned long long)this << std::endl;
		}
		lmbd(ptr);
		initialized = true;
	}

	__host__ void mark_as_initialized(){
		initialized = true;
	}

	template<typename F>
	__host__ void foreach(F&& lmbd){
		for (size_t i = 0; i < N; i++){
			lmbd(ptr[i]);
		}
	}
	
};
/*
#define A_INIT(arr,thisptr,dptr,x) arr.initialize([&](decltype(arr)* thisptr,cpp11::type_identity<decltype(arr)>::type::ptr_type dptr)x); 
#define A_FOREACH(arr,dptr,x) arr.initialize([&](decltype(arr)* thisptr,cpp11::type_identity<decltype(arr)>::type::ptr_type dptr)x); 
*/
template<typename T>
struct CValue :public CArr<T, 1>{
	__device__ __host__ operator T()const{
		return (T)(this->ptr[0]);
	}

	__device__ __host__ CValue& operator=(const T& item){
		this->ptr[0] = item;
		return *this;
	}
	CValue() :CArr<T,1>(){
		dbgprint("CValue default ctor called.\n");

	}
	CValue(T&& init) :CArr<T, 1>(){
		dbgprint("CValue ctor with initial value called.\n");
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


template<typename T,size_t N,size_t Ch,typename Arr=CArr<T,N>>
struct CArrMulti{
private:
	uint8_t c[sizeof(Arr)*Ch];
	devPtr<T> ptr;
public:

	/*
	template<unsigned...I>
	CArrMulti(cpp14::integer_sequence<I...>) :ptr(thrust::device_malloc<T>(N*Ch)), c({Arr(ptr+I*N,true)...}){
		static_assert(std::is_standard_layout<devPtr<T>>::value, " 'device_ptr' of this version of CUDA is not standard layout.\n");
		static_assert(std::is_trivially_copyable<devPtr<T>>::value, " 'device_ptr' of this version of CUDA is not trivially copyable.\n");
		
	}
	*/
	CArrMulti() :ptr(thrust::device_malloc<T>(N*Ch)){

		static_assert(std::is_standard_layout<devPtr<T>>::value, " 'device_ptr' of this version of CUDA is not standard layout.\n");
		static_assert(std::is_trivially_copyable<devPtr<T>>::value, " 'device_ptr' of this version of CUDA is not trivially copyable.\n");
		
		if (!std::is_standard_layout<Arr>::value){
			report_warn([]{
				std::cerr << "Type:" << typeid(Arr).name() << " is not standard layout." << std::endl;
			});
			
		}

		if (!std::is_trivially_copyable<Arr>::value){
			report_warn([]{
				std::cerr << "Type:" << typeid(Arr).name() << " is not trivially copyable." << std::endl;
			});
			
		}
		
		for (size_t i = 0; i < Ch; i++){
			new(reinterpret_cast<void*>(&reinterpret_cast<Arr*>(c)[i])) Arr(ptr + i*N, true);
		}
		dbgprint("CArrMulti ctor called.\n");
	}

	__device__ __host__ Arr& operator[](size_t ch){
		assert(ch < Ch);
		return reinterpret_cast<Arr*>(&c)[ch];
	}

	~CArrMulti(){
		dbgprint("CArrMulti dtor called and freed.\n");
		thrust::device_free(ptr);
	}



};

template<typename T, size_t X, size_t Y, size_t Z, size_t Ch, typename U = CArr3D<T, X, Y, Z>>
struct CArr3DMulti :public CArrMulti<T, X*Y*Z,Ch,U>{
	
};

template<typename T>
using Field3D = CArr3D<T, NX + 1, NY + 1, NZ + 1>;


template<typename T,size_t N>
using Field3DMulti = CArr3DMulti<T, NX + 1, NY + 1, NZ + 1,N,Field3D<T>>; //well-defined Field3D