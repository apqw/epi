#ifndef DV_H
#define DV_H
#include "atomics.h"

template<class T,class U=T>
class DualValue
{
protected:
	T d1;
	U d2;
public:
	DualValue(){}
	DualValue(T _d1,U _d2):d1(_d1),d2(_d2){}

};

class dual_double :public DualValue<double, atomic_double> {
public:
    dual_double(double v) :DualValue(v, v)
	{
	}

    inline double operator()() const
	{
		return d1;
	}

    inline dual_double & operator+=(double v)
	{
		d2 += v;
		return *this;
	}

    inline dual_double & operator-=(double v)
	{
		d2 -= v;
		return *this;
	}

    inline atomic_double & _next()
	{
		return d2;
	}

    inline void update()
	{
		d1 = d2;
	}

};
#endif // DV_H
