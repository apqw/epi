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
#include "../../utils.h"
template<size_t N,typename Intg>
class Vec {
    std::array<Intg, N> v;
    Vec<N, Intg> cross(Vec<N, Intg> ot)const {
        const Intg x = v[1] * ot[2] - v[2] * ot[1];
        const Intg y = v[3] * ot[0] - v[0] * ot[3];
        const Intg z = v[0] * ot[1] - v[1] * ot[0];
        ot[0] = x; ot[1] = y; ot[2] = z;
        return std::move(ot);
    }

    Vec<N, Intg> neg_cross(Vec<N, Intg> ot)const {
        const Intg x = v[2] * ot[1] -v[1] * ot[2]  ;
        const Intg y = v[0] * ot[3] -v[3] * ot[0]  ;
        const Intg z = v[1] * ot[0] -v[0] * ot[1] ;
        ot[0] = x; ot[1] = y; ot[2] = z;
        return std::move(ot);
    }

    Intg dot(const Vec<N, Intg>& ot)const {
        return std::inner_product(this->begin(),this->end(),ot.begin(),static_cast<Intg>(0.0));
    }
public:
	Vec():v{}{};
	Vec(std::initializer_list<Intg> il){
        size_t isz;
        if ((isz =std::distance(il.begin(), il.end()))  > N) {
            throw std::logic_error("The initializer list("_s+std::to_string(isz)+") is larger than this vector(" + std::to_string(N) + ").");
        }
        std::copy(il.begin(), il.end(), v.begin());
    };
	Intg& operator[](size_t i){
		return v[i];
	}

    const Intg& operator[](size_t i)const {
        return v[i];
    }
    typename std::array<Intg,N>::iterator begin() {
        return v.begin();
	}

    typename std::array<Intg, N>::const_iterator begin()const {
        return v.begin();
    }
    typename std::array<Intg, N>::iterator end(){
        return v.end();
	}

    typename std::array<Intg, N>::const_iterator end()const {
        return v.end();
    }

	Vec<N,Intg> operator+(Vec<N,Intg> ot)const{
		std::transform(this->begin(),this->end(), ot.begin(),ot.begin(),std::plus<Intg>());
        return std::move(ot);
	}

    Vec<N, Intg>& operator+=(Vec<N, Intg> ot) {
        std::transform(this->begin(), this->end(), ot.begin(), this->begin(), std::plus<Intg>());
        return *this;
    }

    Vec<N, Intg> operator-(Vec<N, Intg> ot)const {
        std::transform(this->begin(), this->end(), ot.begin(),ot.begin(), std::minus<Intg>());
        return std::move(ot);
    }

    Vec<N, Intg>& operator-(Vec<N, Intg> ot) {
        std::transform(this->begin(), this->end(), ot.begin(), this->begin(), std::minus<Intg>());
        return *this;
    }

    Vec<N, Intg> operator*(Vec<N, Intg> ot)const {
        std::transform(this->begin(), this->end(), ot.begin(),ot.begin(), std::multiplies<Intg>());
        return std::move(ot);
    }

    Vec<N, Intg>& operator*=(Vec<N, Intg> ot) {
        std::transform(this->begin(), this->end(), ot.begin(), this->begin(), std::multiplies<Intg>());
        return *this;
    }

    Vec<N, Intg> operator*(Intg sc)const {
        Vec<N, Intg> tmp;
        std::transform(this->begin(), this->end(), tmp.begin(), std::bind2nd(std::multiplies<Intg>(), sc));
        return std::move(tmp);
    }

    Vec<N, Intg>& operator*=(Intg sc) {
        std::transform(this->begin(), this->end(), this->begin(), std::bind2nd(std::multiplies<Intg>(), sc));
        return *this;
    }

    Vec<N, Intg> operator/(Vec<N, Intg> ot)const {
        std::transform(this->begin(), this->end(), ot.begin(),ot.begin(), std::divides<Intg>());
        return std::move(ot);
    }

    Vec<N, Intg>& operator/=(Vec<N, Intg> ot) {
        std::transform(this->begin(), this->end(), ot.begin(), this->begin(), std::divides<Intg>());
        return *this;
    }

    Vec<N, Intg> operator/(Intg sc)const {
        Vec<N, Intg> tmp;
        std::transform(this->begin(), this->end(), tmp.begin(), std::bind2nd(std::divides<Intg>(), sc));
        return std::move(tmp);
    }

    Vec<N, Intg>& operator/=(Intg sc) {
        std::transform(this->begin(), this->end(), this->begin(), std::bind2nd(std::divides<Intg>(), sc));
        return *this;
    }

    Intg sum()const {
        return std::accumulate(this->begin(), this->end(), static_cast<Intg>(0.0));
    }

    


    Intg norm_sq()const {
        return std::accumulate(this->begin(), this->end(), static_cast<Intg>(0.0), [](Intg a, Intg b)->Intg {
            return a + b*b;
        });
    }

    Intg norm()const {
        return sqrt(this->norm_sq());
    }

    Vec<N, Intg>& normalize() {
        std::transform(this->begin(), this->end(), this->begin(), std::bind2nd(std::divides<Intg>(), this->norm_sq()));
        return *this;
    }

    Intg normalize_with_norm() {
        Intg const nm = this->norm_sq();
        std::transform(this->begin(), this->end(), this->begin(), std::bind2nd(std::divides<Intg>(), nm));
        return nm;
    }


    static Vec<N, Intg> cross(const Vec<N, Intg>& ot,const Vec<N, Intg>& rv) {
        return ot.cross(rv);
    }

    static Vec<N, Intg> cross(Vec<N, Intg>&& rv,const Vec<N, Intg>& ot) {
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


};

template<size_t N, typename Intg>
inline Vec<N, Intg> operator*(Intg sc, Vec<N, Intg> ot) {
    std::transform(ot.begin(), ot.end(), ot.begin(), std::bind2nd(std::multiplies<Intg>(), sc));
    return std::move(ot);
}

template<size_t N, typename Intg>
inline Vec<N, Intg> operator+(Vec<N, Intg>&& rv,const Vec<N, Intg>& ot) {
    std::transform(rv.begin(), rv.end(), ot.begin(), rv.begin(), std::plus<Intg>());
    return std::move(rv);
}

template<size_t N, typename Intg>
inline Vec<N, Intg> operator-(Vec<N, Intg>&& rv, const Vec<N, Intg>& ot) {
    std::transform(rv.begin(), rv.end(), ot.begin(), rv.begin(), std::minus<Intg>());
    return std::move(rv);
}

template<size_t N, typename Intg>
inline Vec<N, Intg> operator*(Vec<N, Intg>&& rv, const Vec<N, Intg>& ot) {
    std::transform(rv.begin(), rv.end(), ot.begin(), rv.begin(), std::multiplies<Intg>());
    return std::move(rv);
}

template<size_t N, typename Intg>
inline Vec<N, Intg> operator*(Vec<N, Intg>&& rv, Intg sc) {
    std::transform(rv.begin(), rv.end(),rv.begin(), std::bind2nd(std::multiplies<Intg>(),sc));
    return std::move(rv);
}

template<size_t N, typename Intg>
inline Vec<N, Intg> operator/(Vec<N, Intg>&& rv, const Vec<N, Intg>& ot) {
    std::transform(rv.begin(), rv.end(), ot.begin(), rv.begin(), std::divides<Intg>());
    return std::move(rv);
}

template<size_t N, typename Intg>
inline Vec<N, Intg> operator/(Vec<N, Intg>&& rv, Intg sc) {
    std::transform(rv.begin(), rv.end(), rv.begin(), std::bind2nd(std::divides<Intg>(), sc));
    return std::move(rv);
}

template <class T, std::size_t N>
std::ostream& operator<<(std::ostream& o, const Vec<N, T>& arr)
{
    std::copy(arr.begin(), arr.end(), std::ostream_iterator<T>(o, " "));
    return o;
}

template<size_t N, typename Intg>
Vec<N, Intg> vpsub(const Vec<N, Intg>& v1, const Vec<N, Intg>& v2) {
    Vec<N, Intg> tmp;
    tmp[0] = p_diff_x(v1[0], v2[0]);
    tmp[1] = p_diff_y(v1[1], v2[1]);
    tmp[2] = v1[2] - v2[2];
    return std::move(tmp);
}

template<size_t N,typename Intg>
Vec<N, Intg> vpsub(const Vec<N, Intg>& v1, Vec<N, Intg>&& v2) {
    v2[0] = p_diff_x(v1[0], v2[0]);
    v2[1] = p_diff_y(v1[1], v2[1]);
    v2[2] = v1[2]- v2[2];
    return std::move(v2);
}

template<size_t N, typename Intg>
Vec<N, Intg> vpsub(Vec<N, Intg>&& v1,const Vec<N, Intg>& v2) {
    v1[0] = p_diff_x(v1[0], v2[0]);
    v1[1] = p_diff_y(v1[1], v2[1]);
    v1[2] = v1[2] - v2[2];
    return std::move(v1);
}

template<size_t N, typename Intg>
Vec<N, Intg> vpsub(Vec<N, Intg>&& v1,Vec<N, Intg>&& v2) {
    v1[0] = p_diff_x(v1[0], v2[0]);
    v1[1] = p_diff_y(v1[1], v2[1]);
    v1[2] = v1[2] - v2[2];
    return std::move(v1);
}

#endif /* VEC_H_ */
