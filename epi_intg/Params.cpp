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

template<>
std::string convert_to_string<bool>(bool t){
	return t?"1" : "0";
}

template<>
std::string convert_to_string<std::string>(std::string t){
	return std::string(str_sign[0])+t+str_sign[0];
}


void Params::s_Ctor() {
    init();
}

void Params::s_Ctor(const std::string& path) {
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
        pi->set_from(kvmap);
    }
}

void Params::generate_paramfile(const std::string& out_path)const {
    std::ofstream of(out_path);
    if (!of) {
        throw std::runtime_error("Failed to create the file:"_s + out_path);
    }

    for (auto& pi : piset) {
        of << pi->name[0] << delim << pi->get_as_string()<<std::endl;
    }


}

void Params::diff(const Params& p1, const Params& p2) {
    size_t pisize = p1.piset.size();
    if (pisize != p2.piset.size()) {
        throw std::logic_error("Comparing different parameter sets");
    }


    std::cout << "Checking parameters..." << std::endl;

    std::cout << "Different parameter(s):\n" << std::endl;
    //std::cout << "<Parameter name>:<Parameter1> <Parameter2>"<< std::endl;
    bool found = false;
    int count = 0;
    for (size_t i = 0; i < pisize; i++) {
        if (!p1.piset[i]->compare_as_data(p2.piset[i])) {
            std::cout << p1.piset[i]->name[0] << ":" << p1.piset[i]->get_as_string() << " " << p2.piset[i]->get_as_string() << std::endl;
            found = true;
            count++;
        }
    }

    
        std::cout <<"\n"<< std::to_string(count)<<" different parameter(s) found." << std::endl;
    

    //std::cout << "Finish." << std::endl;
}

Params::~Params(){
	for(auto& pi:piset){
		delete pi;
	}
}

