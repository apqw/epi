#pragma once
#include "define.h"
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
template<typename TArr,class Func,unsigned int DEPTH>
class __arr_func {
public:
	static void map(TArr& arr) {
		for (auto& a : arr) {
			__arr_func<decltype(a),Func,DEPTH-1>::map(a);
		}
	}
};

template<class Func,typename T>
class __arr_func<T,Func,1> {
public:
	static void map(T& arr) {
		for (auto& a : arr) {
			Func()(a);
		}
	}
};

template<class Func, typename T>
class __arr_func<T, Func, 0> {
public:
	static void map(T& arr) {
//		static_assert(false, "DEPTH should be >0\n");
        assert(!"DEPTH should be >0\n");
	}
};


template<class Func, unsigned int DEPTH, typename TArr>
void arr_map(TArr& arr) {
	__arr_func<TArr, Func, DEPTH>::map(arr);
}

template<class L>
void parallel_for(int rstart,int rend,const L& lmbd) {
	tbb::parallel_for(tbb::blocked_range<int>(rstart, rend), [&](const tbb::blocked_range< int >& range) {
		for (int i = range.begin(); i < range.end(); i++) {
			lmbd(i);
		}
	});
}