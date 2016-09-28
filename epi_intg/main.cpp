/*
 * main.cpp
 *
 *  Created on: 2016/09/20
 *      Author: yasu7890v
 */



#include <iostream>
#include "Params.h"
#include "global.h"
#include "parser.h"
int main(int argc,char** argv){
	Params pa;
	pa.load("/home/yasu7890v/wau/ppsss");
	pa.ANX=1010;

	PW.init_params(pa);

	ParamsWrapper pw;
	pw.init_params(pa);
	//pw.init_params(pa);

	auto pset = parse_paramtext("/home/yasu7890v/wau/ppsss");
	for(auto& kv:pset){
		std::cout<<kv.first<<":"<<kv.second<<std::endl;
	}

	std::cout<<"ok"<<std::endl;
}
