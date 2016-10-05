/*
 * Vec.h
 *
 *  Created on: 2016/10/05
 *      Author: yasu7890v
 */

#ifndef VEC_H_
#define VEC_H_
#include <mkl.h>
#include <initializer_list>
#include <algorithm>
template<size_t N,typename Intg>
class Vec {
	Intg v[N];
	Vec():v{}{};
	Vec(std::initializer_list<Intg> il):v(il){};
	Intg& operator[](size_t i){
		return v[i];
	}
	Intg* begin(){
		return v;
	}
	Intg* end(){
		return &v[N];
	}
	Vec<N,Intg>&& operator+(const Vec<N,Intg>& ot){
		Vec<N,Intg> tmp;
		std::transform(this->begin(),this->end(),ot.begin(),tmp.begin(),std::plus<Intg>());
		return std::move(tmp);
	}
};

#endif /* VEC_H_ */
