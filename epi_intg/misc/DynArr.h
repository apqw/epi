#pragma once
#include <vector>
#include <cstdlib>
#include <type_traits>
#include <cstring>
#include <numeric>
#include <functional>

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

/*

template<size_t Count, typename Last>
int arr_idx_impl(const size_t dim[], Last z) {
    return z;
}

template<size_t Count,typename First,typename Second,typename... Rest>
int arr_idx_impl(const size_t dim[],First f,Second s,Rest... r) {
return f + int(dim[Count]) * arr_idx_impl<Count-1>(dim,s, r...);
}

// forward decl
template<class ...Tn>
struct arr_idx;

// recursion anchor
template<>
struct arr_idx<>
{
    template<size_t Count,class ...Un>
    static int apply(const size_t dim[],Un const&... un)
    {
        return arr_idx_impl<Count>(dim, int(un)...);
    }
};

// recursion
template<class T, class ...Tn>
struct arr_idx<T, Tn...>
{
    template<size_t Count,class ...Un>
    static int apply(const size_t dim[],T const& t, Tn const&... tn, Un const&... un)
    {
        // bubble 1st parameter backwards
        return arr_idx<Tn...>::apply<Count>(dim,tn..., t, un...);
    }
};

template<typename T,typename U,typename Fn>
T accumulate_list(std::initializer_list<T> arr,U init,Fn lmbd) {
    return std::accumulate(arr.begin(), arr.end(), (T)init, lmbd);
}


template<typename T,size_t Dim, typename std::enable_if<!std::is_same<bool, T>::value>::type* = nullptr>
struct DynArr {
private:

    std::vector<T> _data;
public:
    const size_t N[Dim];
    const size_t& NX, NY, NZ, NW;
    size_t layer_size;
    template<typename... Dims,typename std::enable_if<sizeof...(Dims)==Dim>::type* =nullptr>
    DynArr(Dims... _N) :
        N{ _N... }, 
        _data(accumulate_list({(size_t)(_N)...}, 1, std::multiplies<size_t>())),NX(N[0]),NY(N[1]),NZ(N[2]),NW(N[3]) {
        layer_size = _data.size() / N[0];
    }


    template<typename FDim,typename... Dims>
    typename std::vector<T>::reference at(FDim fd,Dims... d) {
        return _data[arr_idx<Dims...>::apply<Dim-1>(N,(int)fd,int(d)...)];
    }

    template<typename FDim, typename... Dims>
    typename std::vector<T>::const_reference at(FDim fd, Dims... d)const {
        return _data[arr_idx<Dims...>::apply<Dim-1>(N, int(fd), int(d)...)];
    }
    void all_memset(int _Val) {
        std::memset(&_data[0], _Val, sizeof(T)*_data.size());
    }
    void x_range_memset(size_t xbegin, size_t xrange, int _Val) {
        std::memset(&_data[0]+xbegin*layer_size, _Val, sizeof(T)*xrange*layer_size);
    }
    void x_layer_memset(size_t z, int _Val) {
        x_range_memset(z, 1, _Val);
    }

};

template<typename T>
using Dyn3DArr = DynArr<T, 3>;

template<typename T>
using Dyn4DArr = DynArr<T, 4>;
*/