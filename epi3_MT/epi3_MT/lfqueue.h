#ifndef LFQUEUE_H
#define LFQUEUE_H
#include <atomic>
#include <array>

template<typename T,unsigned int N>
class LFQueue
{
std::array<T,N> _arr;
std::atomic_uint head;
public:
    LFQueue():head(0){}

    void push_back(T& item){
unsigned int my_head= head++;
_arr[my_head]=item;
    }

    void push_back(const T& item){
unsigned int my_head= head++;
_arr[my_head]=item;
    }
    void clear(){
    head=0;
    }
    const auto& raw_arr() const{
        return _arr;
    }
    unsigned int count()const{
        return head;
    }
};

#endif // LFQUEUE_H
