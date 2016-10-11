#pragma once
#ifndef SWAPDATA_H
#define SWAPDATA_H
#include <unordered_map>
#include <type_traits>
#include <utility>
#include "../utils.h"


template<typename T>
class SwapData
{
    T data1_internal;
    T data2_internal;
    T* data1;
    T* data2;
public:
    SwapData() {
        data1 = &data1_internal;
        data2 = &data2_internal;
    }
    SwapData(T&& data) :data1_internal(data), data2_internal(std::move(data))
    {
        data1 = &data1_internal;
        data2 = &data2_internal;
    }
    SwapData(const T& data) :data1_internal(data), data2_internal(data) {
        data1 = &data1_internal;
        data2 = &data2_internal;
    }

    template<
        typename... U,
        typename std::enable_if<!single_arg_and_eq_ref<T,U...>::value>::type* =0
    >
    SwapData(U... Args) : data1_internal(Args...), data2_internal(Args...) {
        data1 = &data1_internal;
        data2 = &data2_internal;
    }

    T& first() {
        return *data1;
    }
    T& second() {
        return *data2;
    }

    const T& first()const {
        return *data1;
    }
    const T& second()const {
        return *data2;
    }
    void swap() {
        std::swap(data1, data2);
    }
};
#endif