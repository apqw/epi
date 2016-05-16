#ifndef CELL_H
#define CELL_H
#include "define.h"
class CellData
{
private:
    unsigned int cell_num;
    double pos[MAX_CELL_NUM][3];//x,y,z, x,y,z, x,y,z,...
    double velocity[MAX_CELL_NUM][3];
    double radius[MAX_CELL_NUM];
    double age[MAX_CELL_NUM];
    unsigned int neighbor_count[MAX_CELL_NUM];
    unsigned int neighbor_index[MAX_CELL_NUM][MAX_CONNECT_NUM];
    CELL_STATE cell_state[MAX_CELL_NUM];
    bool order_invalidated=false;
public:
    CellData();
    void calc_cell_interaction_force(unsigned int cidx1,unsigned int cidx2,
                                     double* out_fx,double* out_fy,double* out_fz);
    void set_vel_by_cell_interact();
    void set_pos_by_velocity();
    double cell_dist_sq(unsigned int cidx1,unsigned int cidx2);
    void add_cell(CELL_STATE cst,double x=0,double y=0,double z=0,double vx=0,double vy=0,double vz=0,double _rad=rad_default,double _age=0);
    void remove_cell(unsigned int idx);
    void set_neighbor();
};

#endif // CELL_H
