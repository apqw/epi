#pragma once
#include <vector>
#include <cstdlib>
#include <type_traits>
#include <cstring>
template<typename T,typename std::enable_if<!std::is_same<bool,T>::value>::type* =nullptr>
struct Dyn3DArr {
private:

    std::vector<T> _data;
public:
    const size_t NX, NY, NZ;
    Dyn3DArr(size_t _NX, size_t _NY, size_t _NZ) :NX(_NX), NY(_NY), NZ(_NZ),_data(_NX*_NY*_NZ) {}
    typename std::vector<T>::reference at(size_t x, size_t y, size_t z) {
        return _data[z + y*NZ + x*NZ*NY];
    }

    typename std::vector<T>::const_reference at(size_t x, size_t y, size_t z)const {
        return _data[z + y*NZ + x*NZ*NY];
    }
    void all_memset(int _Val) {
        std::memset(&_data[0], _Val, sizeof(T)*NX*NY*NZ);
    }
    void x_range_memset(size_t xbegin,size_t xrange,int _Val){
            	std::memset(&(this->at(xbegin,0,0)), _Val, sizeof(T)*xrange*NY*NZ);
            }
    void x_layer_memset(size_t z,int _Val){
    	x_range_memset(z,1,_Val);
    }


};


template<typename T,typename std::enable_if<!std::is_same<bool,T>::value>::type* =nullptr>
struct Dyn4DArr {
private:

    std::vector<T> _data;
public:
const size_t NX, NY, NZ,NW;
    Dyn4DArr(size_t _NX, size_t _NY, size_t _NZ,size_t _NW) :NX(_NX), NY(_NY), NZ(_NZ),NW(_NW), _data(_NX*_NY*_NZ*_NW) {}
    typename std::vector<T>::reference at(size_t x, size_t y, size_t z,size_t w) {
        return _data[w + z*NW + y*NZ*NW+x*NY*NZ*NW];
    }

    typename std::vector<T>::const_reference at(size_t x, size_t y, size_t z, size_t w)const {
        return _data[w + z*NW + y*NZ*NW+x*NY*NZ*NW];
    }

    void all_memset(int _Val) {
        std::memset(&_data[0], _Val, sizeof(T)*NX*NY*NZ*NW);
    }
};
