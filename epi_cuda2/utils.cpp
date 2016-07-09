/*
 * utilc.cpp
 *
 *  Created on: 2016/07/09
 *      Author: yasu7890v
 */




#include "utils.h"

std::string int_to_string(int number)
{
  std::stringstream ss;
  ss << number;
  return ss.str();
}
