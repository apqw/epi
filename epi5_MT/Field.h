#pragma once
#include <tbb/parallel_for.h>
#include <tbb/blocked_range3d.h>

template<typename T,size_t X,size_t Y,size_t Z>
class Field
{
	T _fdata[X][Y][Z];
	
public:

	Field():_fdata{}{}
	void parallel_body_proc(T&& lmbd) {
		tbb::parallel_for(tbb::blocked_range3d<double>(0,X-1,0,Y-1,0,Z-1), lmbd);
	}

    inline auto&& operator()() {
		return _fdata;
	}
	//~Field();
};

