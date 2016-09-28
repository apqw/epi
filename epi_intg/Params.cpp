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

template<typename T>
T str_cast(const std::string& str){
	return T();
}

template<>
auto str_cast(const std::string& str)->unsigned int{
	return std::stoull(str);
}

Params::Params(){}
void Params::load(const char* path){
	auto kvmap = parse_paramtext(path);
	MEMB_NUM_X =
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

