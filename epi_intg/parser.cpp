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
#include <regex>
#include <array>

static const constexpr std::array<const char*,2> comment_out_signature= {"#","//"};
bool is_commented_out(const std::string& ln){
	for(const char* co:comment_out_signature){
		size_t pos=ln.find(co);
		if(pos==0){
			return true;
		}
	}
	return false;
}
std::string remove_LR_space(const std::string& str) {
    return std::regex_replace(str, std::regex("^ +| +$|( ) +"), "$1");
}
std::map<std::string,std::string> parse_paramtext(const std::string& path){
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

		if(is_commented_out(line))continue;

		if(line.find(delim)==std::string::npos){
			std::cout<<delim<<" not found at line "<<lcount<<"."<<std::endl;
			continue;
		}

		std::istringstream st(line);
		std::getline(st,key,delim);
        key = remove_LR_space(key);
		std::getline(st,rest); //extract the rest ( because this does not include endl)
        rest = remove_LR_space(rest);
		kvset.insert(std::make_pair(key,rest));
	}
	return kvset;
}

