/*
 * define.h
 *
 *  Created on: 2016/09/20
 *      Author: yasu7890v
 */

#ifndef DEFINE_H_
#define DEFINE_H_
#include <cstdint>

#define USE_DOUBLE

#ifdef USE_DOUBLE
typedef double real;
#else
typedef float real;
#endif
enum CELL_STATE:uint32_t {
    ALIVE = 0,
    DEAD = 1,
    DISA = 2,
    UNUSED = 3, //
    FIX = 4,
    BLANK = 5,
    DER = 6,
    MUSUME = 7,
    AIR = 8,
    MEMB = 9
};
static constexpr unsigned int WHOLE = 0;
static constexpr unsigned int BASAL = 1;
#ifdef _WIN32
#define RESTRICT __restrict
#else
#define RESTRICT __restrict__
#endif

#endif /* DEFINE_H_ */
