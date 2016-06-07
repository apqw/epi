#pragma once
#include <tbb/parallel_for.h>
#include <tbb/blocked_range3d.h>
#include <array>

template<typename T,size_t X,size_t Y,size_t Z>
class Field
{
	using Arr = T[X][Y][Z];
	Arr _fdata;
	
public:

	Field():_fdata{}{}
	void parallel_body_proc(T&& lmbd) {
		tbb::parallel_for(tbb::blocked_range3d<double>(0,X-1,0,Y-1,0,Z-1), lmbd);
	}

    inline Arr& operator()() {
		return _fdata;
	}

	inline const Arr& operator()()const {
		return _fdata;
	}
	//~Field();
};

