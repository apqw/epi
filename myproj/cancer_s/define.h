#ifndef DEFINE_H
#define DEFINE_H

static constexpr double LX=50;
static constexpr double LY=50;
static constexpr double LZ=100;
static constexpr int NX=100;
static constexpr int NY=100;
static constexpr int NZ=100;
static constexpr double dx=LX/NX;
static constexpr double dy=LY/NY;
static constexpr double dz=LZ/NZ;


static constexpr int tNX = NX+2;
static constexpr int tNY = NY+2;
static constexpr int tNZ = NZ+2;
static constexpr double dt_n=0.2;

static constexpr double niche_diffuse=0.1;

enum CELL_STATE:unsigned int{
    UNDEFINED=0,
    C_STEM=1u<<1,
    C_PRE=1u<<2,
    C_DIFFED=1u<<3,
    N_STEM=1u<<4,
    N_PRE=1u<<5,
    N_DIFFED=1u<<6
};

#endif // DEFINE_H
