#pragma once
#define cdefd static constexpr double
#define cdefi static constexpr int
#define cdefui static constexpr unsigned int

cdefui MAX_CELL_NUM = 30000;
cdefui MAX_CONNECT_CELL_NUM = 400;

enum CELL_STATE :unsigned int {
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