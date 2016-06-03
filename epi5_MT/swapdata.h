#ifndef SWAPDATA_H
#define SWAPDATA_H

template<typename T>
class SwapData
{
    T data1_internal;
    T data2_internal;
    T* data1;
    T* data2;
public:
    SwapData(){
        data1=&data1_internal;
        data2=&data2_internal;
    }
    SwapData(T&& data):data1_internal(data),data2_internal(data),SwapData(){}
    SwapData(const T& data):data1_internal(data),data2_internal(data),SwapData(){}

    T& first(){
return *data1;
    }
    T& second(){
        return *data2;
    }
    void swap(){
        T* tmp=data2;
        data2=data1;
        data1=tmp;
    }
};

#endif // SWAPDATA_H
