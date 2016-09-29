/*
 * parser.h
 *
 *  Created on: 2016/09/23
 *      Author: yasu7890v
 */

#ifndef PARSER_H_
#define PARSER_H_

#include <map>
#include <string>
static const constexpr char delim = '=';
std::map<std::string,std::string> parse_paramtext(const std::string&);



#endif /* PARSER_H_ */
