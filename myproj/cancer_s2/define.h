#ifndef DEFINE_H
#define DEFINE_H

#define cdefd static constexpr double
#define cdefi static constexpr int
#define cdefui static constexpr unsigned int
#define cceild(x) (((x)<0||(double)((int)(x)) == (x))?((int)(x)):((int)(x)+1))
cdefi MAX_CELL_NUM = 30000;
cdefi MAX_CONNECT_NUM = 400;
enum CELL_STATE{
    UNDEF=0,
    C_STEM=1,
    C_PRE=2,
    C_DIFF=3,
    N_STEM=4,
    N_PRE=5,
    N_DIFF=6
};
cdefd LX = 50.;
cdefd LY = 50.;
cdefd LZ = 100.;

cdefui NX=100;
cdefui NY=100;
cdefui NZ=200;

cdefd dx=LX/NX;
cdefd dy=LY/NY;
cdefd dz=LZ/NZ;

cdefd dt=0.01;
cdefd rad_default=1.4;

cdefd interact_cut_off_range = 5.;

cdefui neighbor_div_x = (int)cceild(LX/interact_cut_off_range);
cdefui neighbor_div_y = (int)cceild(LY/interact_cut_off_range);
cdefui neighbor_div_z = (int)cceild(LZ/interact_cut_off_range);

cdefd neighbor_dx = (double)(LX/neighbor_div_x);
cdefd neighbor_dy = (double)(LY/neighbor_div_y);
cdefd neighbor_dz = (double)(LZ/neighbor_div_z);



cdefd LJ_eps = 0.01;

#endif // DEFINE_H
