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
static const constexpr std::array<const char*, 2> str_sign = { "\"", "\'" };

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
	} catch (std::invalid_argument& ie) {
		throw
		std::invalid_argument("Failed to read the value of "_s + first + "(='" + item + "') with the following exception:\n" + typeid(ie).name() + ": " + ie.what() + "\n");
	} catch (std::out_of_range& oe) {
		throw
		std::out_of_range("Failed to read the value of "_s + first + "(='" + item + "') with the following exception:\n" + typeid(oe).name() + ": " + oe.what() + "\n");
	} catch (std::exception& e) {
		throw
		std::logic_error("Failed to read the value of "_s + first + "(='" + item + "') with the following exception:\n" + typeid(e).name() + ": " + e.what() + "\n");
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
protected:

public:
	std::vector<const char*> name;

	template<typename ... Nm, class = typename std::enable_if<sizeof...(Nm)>=1 && is_same_multiple<const char*,Nm...>::value>::type>
	ParamInfo(Nm... _n):name { {_n...}} {}

	virtual void set_from(const std::map<std::string,std::string>& mp) {}

	virtual std::string get_as_string()const {}
	virtual bool compare_as_data(const ParamInfo* p2) {}
};

template<typename T>
std::string convert_to_string(T t) {
	return std::to_string(t);
}

template<>
std::string convert_to_string<bool>(bool t);

template<>
std::string convert_to_string<std::string>(std::string t);

template<typename T>
struct ParamInfoGeneral: ParamInfo {
private:

	T* __param;
	bool optional;
	typedef ParamInfo Base;
public:
	template<typename ... Nm, class = typename std::enable_if<sizeof...(Nm)>=1 && is_same_multiple<const char*,Nm...>::value>::type>
	ParamInfoGeneral(T& _p, Nm... _n) :__param(&_p),ParamInfo(_n...),optional(false) {}

	template< typename... Nm,class = typename std::enable_if<sizeof...(Nm)>=1 && is_same_multiple<const char*,Nm...>::value>::type>
	ParamInfoGeneral(T& _p, bool _optional,Nm... _n) :ParamInfoGeneral(_p, _n... ) {
		optional=_optional;

	}
	void set_from(const std::map<std::string,std::string>& mp) {
		bool chk = true;
		try {
			*__param=force_vget<T>(mp, name);
		} catch(std::exception& e) {
			if(!optional) {
				throw std::runtime_error("Failed to set parameter with the following reason:\n"_s+e.what());
			}
		}
	}

	std::string get_as_string()const {
		return convert_to_string(*__param);
	}

	bool compare_as_data(const ParamInfo* p2) {
		const ParamInfoGeneral<T>* derivedPtr=dynamic_cast<const ParamInfoGeneral<T>*>(p2);
		if(derivedPtr==nullptr) {
			throw std::logic_error("Comparing different data type");
		}

		return *__param==*derivedPtr->__param;
	}

};

template<typename T, typename ...U>
ParamInfoGeneral<T>* make_paraminfo(T& t, U ...v) {
	return new ParamInfoGeneral<T>(t, v...);
}
#define gpa(name,...) make_paraminfo(name,#name,__VA_ARGS__)
#define gp1(name) make_paraminfo(name,#name)
#define gpo(name) make_paraminfo(name,true,#name)
#define gpoa(name,...) make_paraminfo(name,true,#name,__VA_ARGS__)

class Params {
protected:
	std::vector<ParamInfo*> piset;
	void s_Ctor();
	void s_Ctor(const std::string& paramfile);
public:
	Params();
	Params(const std::string& paramfile);
	virtual void init();
	virtual void generate_paramfile(const std::string& out_path) const;
	virtual void load(const std::string& paramfile);
	static void diff(const Params& p1, const Params& p2);
	virtual ~Params();
};

#endif /* PARAMS_H_ */
