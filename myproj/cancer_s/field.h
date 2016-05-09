#pragma once
#include "define.h"
#include <array>
#include <vector>

template<typename T,unsigned X,unsigned Y,unsigned Z>
using RawArr3D=std::array<std::array<std::array<double,Z>,Y>,X>;
template<typename T>
using Arr3D=std::vector<std::vector<std::vector<T>>>;
template<typename T>
using Arr2D=std::vector<std::vector<T>>;
template<typename T>
using Arr1D=std::vector<T>;



class Field
{
private:
    RawArr3D<double,tNX,tNY,tNZ> map_tmp;
public:
    RawArr3D<double,tNX,tNY,tNZ> niche_value;
    RawArr3D<double,tNX,tNY,tNZ> niche_gen;
    RawArr3D<double,tNX,tNY,tNZ> niche_gen_weight;
    Field():map_tmp{0},niche_value{0},niche_gen{0},niche_gen_weight{0}{}
    void test();
void calc_niche_diffusion();
    template<class L>
    void define_niche_gen(L&& lmbd){
        lmbd(niche_gen,niche_gen_weight);
    }


};

