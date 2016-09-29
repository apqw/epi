#ifndef DV_H
#define DV_H
#include "atomics.h"
#include "../../define.h"
template<class T, class U = T>
class DualValue
{
protected:
    T d1;
    U d2;
public:
    DualValue() {}
    DualValue(T _d1, U _d2) :d1(_d1), d2(_d2) {}

};
class dual_real :public DualValue<real, atomic_fp<real>> {
public:
    dual_real(real v) :DualValue<real, atomic_fp<real>>(v, v)
    {
    }

    inline real operator()() const
    {
        return d1;
    }

    inline dual_real & operator+=(real v)
    {
        d2 += v;
        return *this;
    }

    inline dual_real & operator-=(real v)
    {
        d2 -= v;
        return *this;
    }

    inline atomic_fp<real> & _next()
    {
        return d2;
    }

    inline void update()
    {
        d1 = d2;
    }

    inline void _set(real v) {
        d1 = v;
        d2 = v;
    }

};
#endif // DV_H
