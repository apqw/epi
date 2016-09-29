#pragma once
#include <vector>
#include <cstdlib>
template<typename T>
struct Dyn3DArr {
private:
    size_t NX, NY, NZ;
    std::vector<T> _data;
public:

    Dyn3DArr(size_t _NX, size_t _NY, size_t _NZ) :NX(_NX), NY(_NY), NZ(_NZ),_data(_NX*_NY*_NZ) {}
    T& at(size_t x, size_t y, size_t z) {
        return _data[x + y*NX + z*NX*NY];
    }

    const T& at(size_t x, size_t y, size_t z)const {
        return _data[x + y*NX + z*NX*NY];
    }
    void all_memset(int _Val) {
        std::memset(&_data[0], _Val, sizeof(T)*NX*NY*NZ);
    }
};

template<typename T>
struct Dyn4DArr {
private:
    size_t NX, NY, NZ,NW;
    std::vector<T> _data;
public:

    Dyn4DArr(size_t _NX, size_t _NY, size_t _NZ,size_t _NW) :NX(_NX), NY(_NY), NZ(_NZ),NW(_NW), _data(_NX*_NY*_NZ*_NW) {}
    T& at(size_t x, size_t y, size_t z,size_t w) {
        return _data[x + y*NX + z*NX*NY+w*NX*NY*NZ];
    }

    const T& at(size_t x, size_t y, size_t z, size_t w)const {
        return _data[x + y*NX + z*NX*NY + w*NX*NY*NZ];
    }

    void all_memset(int _Val) {
        std::memset(&_data[0], _Val, sizeof(T)*NX*NY*NZ*NW);
    }
};