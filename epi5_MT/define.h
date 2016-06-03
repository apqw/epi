#ifndef DEFINE_H
#define DEFINE_H

#define cdefd static constexpr double
#define cdefi static constexpr int
#define cdefui static constexpr unsigned int
#define cdefs static constexpr size_t

enum CELL_STATE:size_t {
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

namespace cont{
cdefs MAX_CONNECTED_CELL_NUM=400;
cdefs MAX_CELL_NUM=30000;
cdefd R_max = 1.4;//ok
cdefd R_der = 1.4;//ok
cdefd R_memb = 1.0;//ok
cdefd ca2p_init=0.122;
cdefd IP3_init=0;
cdefd ex_inert_init=0.97;
}

#endif // DEFINE_H
