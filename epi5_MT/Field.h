#pragma once
#include <tbb/parallel_for.h>
#include <tbb/blocked_range3d.h>

template<size_t X,size_t Y,size_t Z>
class Field
{
	double _fdata[X+1][Y+1][Z+1];
	
public:

	Field():_fdata{}{}
	template<typename T>
	void parallel_body_proc(T&& lmbd) {
		tbb::parallel_for(tbb::blocked_range3d<double>(0,X,0,Y,0,Z), lmbd);
	}

	//~Field();
};

