/*
 * parser.cpp
 *
 *  Created on: 2016/09/23
 *      Author: yasu7890v
 */


#include "parser.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cstdlib>
static const constexpr char delim='=';
std::map<std::string,std::string> parse_paramtext(const char* path){
	std::ifstream tf(path);
	if(!tf){
		std::cerr<<"Parse failed."<<path<<std::endl;
		std::exit(1);
	}
	std::string line;
	std::string key,rest;
	std::map<std::string,std::string> kvset;
	int lcount=0;
	while(std::getline(tf,line)){
		lcount++;
		if(line.find(delim)==std::string::npos){
			std::cout<<delim<<" not found at line "<<lcount<<"."<<std::endl;
			continue;
		}

		std::istringstream st(line);
		std::getline(st,key,delim);
		std::getline(st,rest); //extract the rest ( because this does not include endl)
		kvset.insert(std::make_pair(key,rest));
	}
	return kvset;
}

