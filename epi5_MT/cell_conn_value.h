#ifndef CELL_CONN_VALUE_H
#define CELL_CONN_VALUE_H
#include "define.h"
template<typename T>
class cell_conn_value
{
protected:
    T _value;
    bool _conn_checked=false;
public:
    cell_conn_value(){}
    cell_conn_value(T&& value):_value(value){}
    cell_conn_value(const T& value):_value(value){}
    inline void check(){
        _conn_checked=true;
    }
    inline void uncheck(){
        _conn_checked=false;
    }

    inline bool is_checked(){
        return _conn_checked;
    }

    inline T&& operator()(){
        return _value;
    }
};

class gj_value:public cell_conn_value<double>{
public:
    gj_value():cell_conn_value(cont::gj_init){}
    gj_value(double d):cell_conn_value(d){}
    inline void reset(){
        _value=cont::gj_init;
    }
};

#endif // CELL_CONN_VALUE_H
