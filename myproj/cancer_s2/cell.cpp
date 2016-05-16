#include "cell.h"
#include "funcs.h"
#include <cassert>
#include <cstring>
CellData::CellData()
{

}

double CellData::cell_dist_sq(unsigned int cidx1,unsigned int cidx2){
return p_dist_sq(pos[cidx1][0],pos[cidx1][1],pos[cidx1][2],pos[cidx2][0],pos[cidx2][1],pos[cidx2][2]);
}

void CellData::calc_cell_interaction_force(unsigned int cidx1, unsigned int cidx2,
                                    double *out_fx, double *out_fy, double *out_fz){
double diffx=p_diff_sc_x(pos[cidx1][0],pos[cidx2][0]);
double diffy=p_diff_sc_y(pos[cidx1][1],pos[cidx2][1]);
double diffz=pos[cidx1][2]-pos[cidx2][2];
double fc = calc_LJ_force_coef(diffx*diffx+diffy*diffy+diffz*diffz,radius[cidx1],radius[cidx2]);
*out_fx=fc*diffx;
*out_fy=fc*diffy;
*out_fz=fc*diffz;
}

void CellData::set_vel_by_cell_interact(){

for(int i=0;i<cell_num;++i){
    double fx,fy,fz;
    for(int j=0;j<i;++j){

        calc_cell_interaction_force(i,j,&fx,&fy,&fz);
        velocity[i][0]+=dt*fx;
        velocity[i][1]+=dt*fy;
        velocity[i][2]+=dt*fz;

        velocity[j][0]-=dt*fx;
        velocity[j][1]-=dt*fy;
        velocity[j][2]-=dt*fz;
    }

}
}

void CellData::set_pos_by_velocity(){
for(int i=0;i<cell_num;++i){
    pos[i][0]+=dt*velocity[i][0];
    pos[i][1]+=dt*velocity[i][1];
    pos[i][2]+=dt*velocity[i][2];
}
}

void CellData::add_cell(CELL_STATE cst,double x=0,double y=0,double z=0,double vx=0,double vy=0,double vz=0,double _rad=rad_default,double _age=0){
    pos[cell_num][0]=x;
    pos[cell_num][1]=y;
    pos[cell_num][2]=z;

    velocity[cell_num][0]=vx;
    velocity[cell_num][1]=vy;
    velocity[cell_num][2]=vz;

    radius[cell_num]=_rad;
    age[cell_num]=_age;

    cell_state[cell_num]=cst;
    neighbor_count[cell_num]=0;
    ++cell_num;
}

void CellData::remove_cell(unsigned int idx){
    assert(idx<cell_num&&cell_num>0);
    order_invalidated=true;
    --cell_num;
    pos[idx][0]=pos[cell_num][0];
    pos[idx][1]=pos[cell_num][1];
    pos[idx][2]=pos[cell_num][2];

    velocity[idx][0]=velocity[cell_num][0];
    velocity[idx][1]=velocity[cell_num][1];
    velocity[idx][2]=velocity[cell_num][2];

    radius[idx]=radius[cell_num];
    age[idx]=age[cell_num];
    cell_state[idx]=cell_state[cell_num];
    neighbor_count[idx]=neighbor_count[cell_num];
    std::memcpy(&(neighbor_index[idx]),&(neighbor_index[cell_num]),sizeof(unsigned int)*MAX_CONNECT_NUM);
}

void CellData::set_neighbor(){
static unsigned int neighbor_cell_index[neighbor_div_x][neighbor_div_y][neighbor_div_z][MAX_CONNECT_NUM];

}
