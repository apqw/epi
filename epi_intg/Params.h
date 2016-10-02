/*
 * Params.h
 *
 *  Created on: 2016/09/20
 *      Author: yasu7890v
 */
#pragma once
#ifndef PARAMS_H_
#define PARAMS_H_
#include "define.h"
#include "utils.h"
#include <vector>
#include <type_traits>
#include <typeinfo>
#include <initializer_list>
#include <map>
#include <stdexcept>
#include <string>
#include <array>
static const constexpr std::array<const char*, 2> str_sign = { "\"","\'" };

template<typename T>
T str_cast(const std::string& str) {
    return T();
}

template<>
unsigned int str_cast<unsigned int>(const std::string& str);

template<>
int str_cast<int>(const std::string& str);

template<>
double str_cast<double>(const std::string& str);

template<>
float str_cast<float>(const std::string& str);

template<>
bool str_cast<bool>(const std::string& str);

template<>
std::string str_cast<std::string>(const std::string& str);

template<typename VarTy>
VarTy try_get_with_cast(const std::string& item, const char* first = "NONE") {
    try {
        return str_cast<VarTy>(item);
    }
    catch (std::invalid_argument& ie) {
        throw std::invalid_argument("Failed to read the value of "_s + first + "(='" + item + "') with the following exception:\n" + typeid(ie).name() + ": " + ie.what() + "\n");
    }
    catch (std::out_of_range& oe) {
        throw std::out_of_range("Failed to read the value of "_s + first + "(='" + item + "') with the following exception:\n" + typeid(oe).name() + ": " + oe.what() + "\n");
    }
    catch (std::exception& e) {
        throw std::logic_error("Failed to read the value of "_s + first + "(='" + item + "') with the following exception:\n" + typeid(e).name() + ": " + e.what() + "\n");
    }
}

template<typename VarTy, typename K, typename V, typename K2>
VarTy force_vget(const std::map<K, V>& mp, const std::vector<K2>& kv) {
    for (auto& key : kv) {
        auto it = mp.find(key);
        if (it != mp.end()) {
            return try_get_with_cast<VarTy>(it->second, key);
        }
    }
    throw std::logic_error(std::string(kv[0]) + " must be specified.");
}

struct ParamInfo {
private:
    enum TypeEnum {
        UINT, INT, REAL,BOOL,STRING
    };


    void* __param;
    TypeEnum pTy;
    template<typename Ty, typename Nm>
    ParamInfo(Ty& _p, std::initializer_list<Nm>_n, typename std::enable_if<std::is_same<Ty, unsigned int>::value>::type* = nullptr) :__param(&_p), name(_n), pTy(TypeEnum::UINT) {
    }

    template<typename Ty, typename Nm>
    ParamInfo(Ty& _p, std::initializer_list<Nm>_n, typename std::enable_if<std::is_same<Ty, int>::value>::type* = nullptr) : __param(&_p), name(_n), pTy(TypeEnum::INT) {
    }
    template<typename Ty, typename Nm>
    ParamInfo(Ty& _p, std::initializer_list<Nm>_n, typename std::enable_if<std::is_same<Ty, real>::value>::type* = nullptr) : __param(&_p), name(_n), pTy(TypeEnum::REAL) {
    }
    template<typename Ty, typename Nm>
    ParamInfo(Ty& _p, std::initializer_list<Nm>_n, typename std::enable_if<std::is_same<Ty, bool>::value>::type* = nullptr) : __param(&_p), name(_n), pTy(TypeEnum::BOOL) {
    }

    template<typename Ty, typename Nm>
    ParamInfo(Ty& _p, std::initializer_list<Nm>_n, typename std::enable_if<std::is_same<Ty, std::string>::value>::type* = nullptr) : __param(&_p), name(_n), pTy(TypeEnum::STRING) {
    }

public:
    std::vector<const char*> name;

    
    
    template<typename Ty, typename... Nm,class = typename std::enable_if<sizeof...(Nm)>=1>::type>
    ParamInfo(Ty& _p, Nm... _n) :ParamInfo(_p, { _n... }) {}

    template<typename K, typename V>
    void set_from(const std::map<K, V>& mp) {
        bool chk = true;
        switch (pTy) {
        case UINT:

            *reinterpret_cast<unsigned int*>(__param) = force_vget<unsigned int>(mp, name);
            break;
        case INT:
            *reinterpret_cast<int*>(__param) = force_vget<int>(mp, name);
            break;
        case REAL:
            *reinterpret_cast<real*>(__param) = force_vget<real>(mp, name);
            break;
        case BOOL:
            *reinterpret_cast<bool*>(__param) = force_vget<bool>(mp, name);
            break;
        case STRING:
            *reinterpret_cast<std::string*>(__param) = force_vget<std::string>(mp, name);
            break;
        default:
            break;
        }
    }

    std::string get_as_string()const {
        bool chk = true;
        switch (pTy) {
        case UINT:
            return std::to_string(*reinterpret_cast<unsigned int*>(__param));
            break;
        case INT:
            return std::to_string(*reinterpret_cast<int*>(__param));
            break;
        case REAL:
            return std::to_string(*reinterpret_cast<real*>(__param));
            break;
        case BOOL:
            return *reinterpret_cast<bool*>(__param) ? "1" : "0";
            break;
        case STRING:
            return std::string(str_sign[0])+*reinterpret_cast<std::string*>(__param)+str_sign[0];
            break;
        default:
            return "!ERROR";
            break;
        }
    }

    static bool compare_as_data(const ParamInfo& p1, const ParamInfo& p2);

};

#define gpa(name,...) ParamInfo(name,#name,__VA_ARGS__)
#define gp1(name) ParamInfo(name,#name)

class Params {
protected:
    std::vector<ParamInfo> piset;
    void s_Ctor();
    void s_Ctor(const std::string& paramfile);
public:
	 Params();
     Params(const std::string& paramfile);
     virtual void init();
     virtual void generate_paramfile(const char* out_path)const;
	 virtual void load(const std::string& paramfile);
     static void diff(const Params& p1, const Params& p2);
};



#endif /* PARAMS_H_ */
