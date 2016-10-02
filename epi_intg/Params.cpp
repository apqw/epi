/*
 * Params.cpp
 *
 *  Created on: 2016/09/20
 *      Author: yasu7890v
 */

#include "Params.h"
#include "parser.h"
#include <iostream>

#include <fstream>


#include <cassert>
#include <regex>


template<>
unsigned int str_cast<unsigned int>(const std::string& str) {
    return std::stoul(str);
}

template<>
int str_cast<int>(const std::string& str) {
    return std::stoi(str);
}

template<>
double str_cast<double>(const std::string& str) {
    return std::stod(str);
}

template<>
float str_cast<float>(const std::string& str) {
    return std::stof(str);
}

template<>
bool str_cast<bool>(const std::string& str) {
    return std::stoi(str) == 1;
}

template<>
std::string str_cast<std::string>(const std::string& str) {

    for (const char* s : str_sign) {
        std::regex re("^"_s + s + "(.*)" + s + "$");
        std::smatch match;
        if (std::regex_match(str, match, re)) {
            return match[1];
        }
    }
    throw std::logic_error("The following string data:" + str + " is ill-formatted.");
}


void Params::s_Ctor() {
    init();
}

void Params::s_Ctor(const std::string& path) {
    std::cout << "super s_ctor w path" << std::endl;
    init();
    const_cast<Params*>(this)->load(path);
}



Params::Params() {
    //s_Ctor();
}

Params::Params(const std::string& path) {
    //s_Ctor(path);
}

void Params::init() {
    throw std::logic_error("init() is not implemented.");
}
void Params::load(const std::string& path){
	auto kvmap = parse_paramtext(path);
    for (auto& pi : piset) {
        pi.set_from(kvmap);
    }
}

void Params::generate_paramfile(const char* out_path)const {
    std::ofstream of(out_path);
    if (!of) {
        throw std::runtime_error("Failed to create the file:"_s + out_path);
    }

    for (auto& pi : piset) {
        of << pi.name[0] << delim << pi.get_as_string()<<std::endl;
    }


}


