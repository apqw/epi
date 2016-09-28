/*
 * Params.cpp
 *
 *  Created on: 2016/09/20
 *      Author: yasu7890v
 */

#include "Params.h"
#include "parser.h"

#include <iostream>
#include <cstdlib>
#include <stdexcept>

#include <type_traits>
#include <typeinfo>
#include <initializer_list>
#include <vector>
#include <cassert>
template<typename T>
T str_cast(const std::string& str){
	return T();
}

template<>
unsigned int str_cast<unsigned int>(const std::string& str){
	return std::stoul(str);
}

template<>
int str_cast<int>(const std::string& str){
	return std::stoi(str);
}

template<>
double str_cast<double>(const std::string& str){
	return std::stod(str);
}

template<>
float str_cast<float>(const std::string& str){
	return std::stof(str);
}

template<typename VarTy>
VarTy try_get_with_cast(const std::string& item,const char* first="NONE"){
	try{
	return str_cast<VarTy>(item);
	}catch(std::invalid_argument& ie){
		throw std::invalid_argument(std::string("Failed to read the value of ")+first+"(='"+item+"') with the following exception:\n"+typeid(ie).name()+": "+ie.what()+"\n");
	}catch(std::out_of_range& oe){
		throw std::out_of_range(std::string("Failed to read the value of ")+first+"(='"+item+"') with the following exception:\n"+typeid(oe).name()+": "+oe.what()+"\n");
	}catch(std::exception& e){
		throw std::logic_error(std::string("Failed to read the value of ")+first+"(='"+item+"') with the following exception:\n"+typeid(e).name()+": "+e.what()+"\n");
	}
}

template<typename VarTy,typename K,typename V,typename K2>
VarTy force_vget(const std::map<K,V>& mp,const std::vector<K2>& kv){
	for(auto& key:kv){
		auto it=mp.find(key);
		if(it!=mp.end()){
			return try_get_with_cast<VarTy>(it->second,key);
		}
	}
	throw std::logic_error(std::string(kv[0])+" must be specified.");
}


struct ParamInfo{
private:
enum TypeEnum{
	UINT,INT,REAL
};

	void* __param;
	TypeEnum pTy;
	template<typename Ty,typename Nm>
	ParamInfo(Ty& _p,std::initializer_list<Nm>_n,typename std::enable_if<std::is_same<Ty,unsigned int>::value>::type* =nullptr):__param(&_p),name(_n),pTy(TypeEnum::UINT){
	}

	template<typename Ty,typename Nm>
		ParamInfo(Ty& _p,std::initializer_list<Nm>_n,typename std::enable_if<std::is_same<Ty,int>::value>::type* =nullptr):__param(&_p),name(_n),pTy(TypeEnum::INT){
		}
	template<typename Ty,typename Nm>
		ParamInfo(Ty& _p,std::initializer_list<Nm>_n,typename std::enable_if<std::is_same<Ty,real>::value>::type* =nullptr):__param(&_p),name(_n),pTy(TypeEnum::REAL){
		}
public:
	std::vector<const char*> name;
	template<typename Ty,typename... Nm>
ParamInfo(Ty& _p,Nm... _n):ParamInfo(_p,{_n...}){}



	template<typename K,typename V>
	void set_from(const std::map<K,V>& mp){
		bool chk=true;
		switch(pTy){
		case UINT:

			*reinterpret_cast<unsigned int*>(__param)=force_vget<unsigned int>(mp,name);
			break;
		case INT:
			*reinterpret_cast<int*>(__param)=force_vget<int>(mp,name);
						break;
		case REAL:
						*reinterpret_cast<real*>(__param)=force_vget<real>(mp,name);
						break;
		default:
			break;
		}
	}

};
#define gpa(name,...) ParamInfo(name,#name,__VA_ARGS__)
#define gp1(name) ParamInfo(name,#name)

Params::Params(){}
void Params::load(const char* path){
	static std::vector<ParamInfo> piset={
		gpa(MEMB_NUM_X,"NMX"),
		gpa(MEMB_NUM_Y,"NMY"),
		gp1(MAX_CONNECT_CELL_NUM),
		gp1(DT_Cell),
		gp1(DT_Ca),
		gpa(MALIG_NUM,"MALIGNANT"),
		gp1(LX),gp1(LY),gp1(LZ),
		gp1(NX),gp1(NY),gp1(NZ),
		gp1(R_max),gp1(R_der),gp1(R_memb),
		gp1(COMPRESS_FACTOR),
		gp1(MEMB_ADHE_RANGE)
	};
	auto kvmap = parse_paramtext(path);

for(auto& pi:piset){
	pi.set_from(kvmap);
}
}

ParamsWrapper::ParamsWrapper():param_is_set(false){}
void ParamsWrapper::init_params(const Params& p){
	std::cout<<"Initializing parameters..."<<std::endl;
	if(param_is_set){
		std::cerr<<"An initialization of the parameters is attempted but the parameters are already set."<<std::endl;
		std::exit(1);
	}
	std::cout<<"Done."<<std::endl;
	*const_cast<Params*>(&params)=p;
	param_is_set=true;
}

