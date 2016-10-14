/*
* Vec.h
*
*  Created on: 2016/10/05
*      Author: yasu7890v
*/

#ifndef VEC_H_
#define VEC_H_
#include <initializer_list>
#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <numeric>
#include <ostream>
#include <iterator>
#include <iostream>
#include "../../utils.h"




template<size_t N, typename Intg>
class Vec {
    std::array<Intg, N> v;
    Vec<N, Intg> cross(Vec<N, Intg> ot)const {
        const Intg x = v[1] * ot.v[2] - v[2] * ot.v[1];
        const Intg y = v[2] * ot.v[0] - v[0] * ot.v[2];
        const Intg z = v[0] * ot.v[1] - v[1] * ot.v[0];
        ot.v[0] = x; ot.v[1] = y; ot.v[2] = z;
        return std::move(ot);
    }

    Vec<N, Intg> neg_cross(Vec<N, Intg> ot)const {
        const Intg x = v[2] * ot.v[1] - v[1] * ot.v[2];
        const Intg y = v[0] * ot.v[2] - v[2] * ot.v[0];
        const Intg z = v[1] * ot.v[0] - v[0] * ot.v[1];
        ot.v[0] = x; ot.v[1] = y; ot.v[2] = z;
        return std::move(ot);
    }

    Intg dot(const Vec<N, Intg>& ot)const {
        Intg t = static_cast<Intg>(0.0);
        for (size_t i = 0; i < N; i++) {
            t += v[i] * ot.v[i];
        }
        return t;
        //return std::inner_product(this->begin(),this->end(),ot.begin(),static_cast<Intg>(0.0));
    }
public:
    Vec() noexcept:v{} {};

    template<typename...U>
    Vec(U...elem) noexcept:v{ elem... } {
        /*
        size_t isz;
        if ((isz =std::distance(il.begin(), il.end()))  > N) {
        throw std::logic_error("The initializer list("_s+std::to_string(isz)+") is larger than this vector(" + std::to_string(N) + ").");
        }
        */
        //std::copy(il.begin(), il.end(), v.begin());

        /*
        for (int i = 0; i < N; i++) {
        v[i] = *(il.begin() + i);
        }
        */
    };
    Intg& operator[](size_t i) {
        return v[i];
    }

    const Intg& operator[](size_t i)const {
        return v[i];
    }
    typename std::array<Intg, N>::iterator begin() {
        return v.begin();
    }

    typename std::array<Intg, N>::const_iterator begin()const {
        return v.begin();
    }
    typename std::array<Intg, N>::iterator end() {
        return v.end();
    }

    typename std::array<Intg, N>::const_iterator end()const {
        return v.end();
    }

    Vec<N, Intg> operator+(Vec<N, Intg> ot)const {
        for (size_t i = 0; i < N; i++) {
            ot.v[i] += v[i];
        }
        return std::move(ot);
    }

    Vec<N, Intg>& operator+=(Vec<N, Intg> ot) {
        for (size_t i = 0; i < N; i++) {
            v[i] += ot.v[i];
        }
        return *this;
    }

    Vec<N, Intg> operator-(Vec<N, Intg> ot)const {
        for (size_t i = 0; i < N; i++) {
            ot.v[i] = v[i] - ot.v[i];
        }
        return std::move(ot);
    }

    Vec<N, Intg>& operator-(Vec<N, Intg> ot) {
        for (size_t i = 0; i < N; i++) {
            v[i] -= ot.v[i];
        }
        return *this;
    }

    Vec<N, Intg> operator*(Vec<N, Intg> ot)const {
        for (size_t i = 0; i < N; i++) {
            ot.v[i] *= v[i];
        }
        return std::move(ot);
    }

    Vec<N, Intg>& operator*=(Vec<N, Intg> ot) {
        for (size_t i = 0; i < N; i++) {
            v[i] *= ot.v[i];
        }
        return *this;
    }

    Vec<N, Intg> operator*(Intg sc)const {
        Vec<N, Intg> tmp;
        for (size_t i = 0; i < N; i++) {
            tmp[i] = sc*v[i];
        }
        return std::move(tmp);
    }

    Vec<N, Intg>& operator*=(Intg sc) {
        for (size_t i = 0; i < N; i++) {
            v[i] *= sc;
        }
        return *this;
    }

    Vec<N, Intg> operator/(Vec<N, Intg> ot)const {
        for (size_t i = 0; i < N; i++) {
            ot.v[i] = v[i] / ot.v[i];
        }
        // std::transform(this->begin(), this->end(), ot.begin(),ot.begin(), std::divides<Intg>());
        return std::move(ot);
    }

    Vec<N, Intg>& operator/=(Vec<N, Intg> ot) {
        const Intg nm = this->norm();
        for (size_t i = 0; i < N; i++) {
            v[i] /= ot.v[i];
        }
        return *this;
    }

    Vec<N, Intg> operator/(Intg sc)const {
        Vec<N, Intg> tmp;
        for (size_t i = 0; i < N; i++) {
            tmp[i] = v[i] / sc;
        }
        // std::transform(this->begin(), this->end(), tmp.begin(), std::bind2nd(std::divides<Intg>(), sc));
        return std::move(tmp);
    }

    Vec<N, Intg>& operator/=(Intg sc) {
        for (size_t i = 0; i < N; i++) {
            v[i] /= sc;
        }
        return *this;
    }

    Intg sum()const {
        Intg t = static_cast<Intg>(0.0);
        for (size_t i = 0; i < N; i++) {
            t += v[i];
        }
        return t;
        //return std::accumulate(this->begin(), this->end(), static_cast<Intg>(0.0));
    }




    Intg norm_sq()const {

        Intg t = static_cast<Intg>(0.0);
        for (size_t i = 0; i < N; i++) {
            t += v[i] * v[i];
        }
        return t;

        /*
        return std::accumulate(this->begin(), this->end(), static_cast<Intg>(0.0), [](Intg a, Intg b)->Intg {
        return a + b*b;
        });
        */
    }

    Intg norm()const {
        return sqrt(this->norm_sq());
    }

    Vec<N, Intg>& normalize() {

        const Intg nm = this->norm();
        for (size_t i = 0; i < N; i++) {
            v[i] /= nm;
        }

        // std::transform(this->begin(), this->begin()+N, this->begin(), std::bind2nd(std::divides<Intg>(), this->norm()));
        return *this;
    }

    Intg normalize_with_norm() {
        const Intg nm = this->norm();
        for (size_t i = 0; i < N; i++) {
            v[i] /= nm;
        }
        return nm;
    }


    static Vec<N, Intg> cross(const Vec<N, Intg>& ot, const Vec<N, Intg>& rv) {
        return ot.cross(rv);
    }

    static Vec<N, Intg> cross(Vec<N, Intg>&& rv, const Vec<N, Intg>& ot) {
        return ot.neg_cross(rv);
    }

    static Vec<N, Intg> cross(const Vec<N, Intg>& ot, Vec<N, Intg>&& rv) {
        return ot.cross(rv);
    }

    static Vec<N, Intg> cross(Vec<N, Intg>&& ot, Vec<N, Intg>&& rv) {
        return ot.cross(rv);
    }

    static Intg dot(const Vec<N, Intg>& ot, const Vec<N, Intg>& rv) {
        return ot.dot(rv);
    }

    template<size_t _N,typename _Intg>
    friend Vec<_N, _Intg> operator*(_Intg sc, Vec<_N, _Intg> ot);

    template<size_t _N, typename _Intg>
    friend Vec<_N, _Intg> operator+(Vec<_N, _Intg>&& rv, const Vec<_N, _Intg>& ot);

    template<size_t _N, typename _Intg>
    friend Vec<_N, _Intg> operator-(Vec<_N, _Intg>&& rv, const Vec<_N, _Intg>& ot);

    template<size_t _N, typename _Intg>
    friend Vec<_N, _Intg> operator*(Vec<_N, _Intg>&& rv, const Vec<_N, _Intg>& ot);

    template<size_t _N, typename _Intg>
    friend Vec<_N, _Intg> operator*(Vec<_N, _Intg>&& rv, _Intg sc);

    template<size_t _N, typename _Intg>
    friend Vec<_N, _Intg> operator/(Vec<_N, _Intg>&& rv, const Vec<_N, _Intg>& ot);

    template<size_t _N, typename _Intg>
    friend Vec<_N, _Intg> operator/(Vec<_N, _Intg>&& rv, _Intg sc);

    template <class T, std::size_t _N>
    friend std::ostream& operator<<(std::ostream& o, const Vec<_N, T>& arr);
};

template<size_t N, typename Intg>
inline Vec<N, Intg> operator*(Intg sc, Vec<N, Intg> ot) {
    for (size_t i = 0; i < N; i++) {
        ot.v[i] *= sc;
    }
    //std::transform(ot.begin(), ot.end(), ot.begin(), std::bind2nd(std::multiplies<Intg>(), sc));
    return std::move(ot);
}

template<size_t N, typename Intg>
inline Vec<N, Intg> operator+(Vec<N, Intg>&& rv, const Vec<N, Intg>& ot) {
    for (size_t i = 0; i < N; i++) {
        rv[i] += ot.v[i];
    }
    return std::move(rv);
}

template<size_t N, typename Intg>
inline Vec<N, Intg> operator-(Vec<N, Intg>&& rv, const Vec<N, Intg>& ot) {
    for (size_t i = 0; i < N; i++) {
        rv[i] -= ot.v[i];
    }
    //std::transform(rv.begin(), rv.end(), ot.begin(), rv.begin(), std::minus<Intg>());
    return std::move(rv);
}

template<size_t N, typename Intg>
inline Vec<N, Intg> operator*(Vec<N, Intg>&& rv, const Vec<N, Intg>& ot) {
    for (size_t i = 0; i < N; i++) {
        rv[i] *= ot.v[i];
    }
    return std::move(rv);
}

template<size_t N, typename Intg>
inline Vec<N, Intg> operator*(Vec<N, Intg>&& rv, Intg sc) {
    for (size_t i = 0; i < N; i++) {
        rv[i] *= sc;
    }
    return std::move(rv);
}

template<size_t N, typename Intg>
inline Vec<N, Intg> operator/(Vec<N, Intg>&& rv, const Vec<N, Intg>& ot) {
    for (size_t i = 0; i < N; i++) {
        rv[i] /= ot.v[i];
    }
    return std::move(rv);
}

template<size_t N, typename Intg>
inline Vec<N, Intg> operator/(Vec<N, Intg>&& rv, Intg sc) {
    for (size_t i = 0; i < N; i++) {
        rv[i] /= sc;
    }
    return std::move(rv);
}

template <class T, std::size_t N>
inline std::ostream& operator<<(std::ostream& o, const Vec<N, T>& arr)
{
    std::copy(arr.begin(), arr.end(), std::ostream_iterator<T>(o, " "));
    return o;
}


template<size_t N, typename Intg>
inline Vec<N, Intg> vpsub(const Vec<N, Intg>& v1, const Vec<N, Intg>& v2) {
    return Vec<N, Intg>(
        p_diff_x(v1[0], v2[0]) ,
        p_diff_y(v1[1], v2[1]) ,
        v1[2] - v2[2]
    );
}


template<size_t N, typename Intg>
inline Vec<N, Intg> vpsub(const Vec<N, Intg>& v1, Vec<N, Intg>&& v2) {
    v2[0] = p_diff_x(v1[0], v2[0]);
    v2[1] = p_diff_y(v1[1], v2[1]);
    v2[2] = v1[2] - v2[2];
    return std::move(v2);
}

template<size_t N, typename Intg>
inline Vec<N, Intg> vpsub(Vec<N, Intg>&& v1, const Vec<N, Intg>& v2) {
    v1[0] = p_diff_x(v1[0], v2[0]);
    v1[1] = p_diff_y(v1[1], v2[1]);
    v1[2] = v1[2] - v2[2];
    return std::move(v1);
}

template<size_t N, typename Intg>
inline Vec<N, Intg> vpsub(Vec<N, Intg>&& v1, Vec<N, Intg>&& v2) {
    v1[0] = p_diff_x(v1[0], v2[0]);
    v1[1] = p_diff_y(v1[1], v2[1]);
    v1[2] = v1[2] - v2[2];
    return std::move(v1);
}


#endif /* VEC_H_ */
