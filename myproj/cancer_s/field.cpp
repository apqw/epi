#include "field.h"
#include <cstring>
#include <iostream>
#include <omp.h>
void Field::test(){
    printf("oo%d",map_tmp.size());
}

void Field::calc_niche_diffusion(){

    for(int i=0;i<tNX;i++){
        for(int j=0;j<tNY;j++){
            for(int k=0;k<tNZ;k++){
                //alpha blending
                //weight‚Í0.0‚Æ1.0‚µ‚©Žg‚í‚È‚¢—\’è
                map_tmp[i][j][k]=niche_value[i][j][k];

            }
        }
    }




//0-y-z,tNX-y-z border plane
//perio
for(int i=0;i<tNY;i++){
    for(int j=0;j<tNZ;j++){
        map_tmp[0][i][j]=map_tmp[tNX-2][i][j];
        map_tmp[tNX-1][i][j]=map_tmp[1][i][j];
    }
}

//x-0-z,x-tNY-z border plane
//perio
for(int i=0;i<tNX;i++){
    for(int j=0;j<tNZ;j++){
        map_tmp[i][0][j]=map_tmp[i][tNY-2][j];
        map_tmp[i][tNY-1][j]=map_tmp[i][1][j];
    }
}
//corner
for(int i=0;i<tNZ;i++){
    map_tmp[0][0][i]=map_tmp[tNX-2][tNY-2][i];
    map_tmp[tNX-1][tNY-1][i]=map_tmp[1][1][i];
    map_tmp[0][tNY-1][i]=map_tmp[tNX-2][1][i];
    map_tmp[tNX-1][0][i]=map_tmp[1][tNY-2][i];
}

        //x-y-0,x-y-tNZ border plane
        //neumann
        for(int i=0;i<tNX;i++){
            for(int j=0;j<tNY;j++){
                map_tmp[i][j][0]=map_tmp[i][j][1];
                map_tmp[i][j][tNZ-1]=map_tmp[i][j][tNZ-2];
            }
        }

for(int i=1;i<=NX;i++){
    for(int j=1;j<=NY;j++){
        for(int k=1;k<=NZ;k++){
            niche_value[i][j][k]+=dt_n*niche_diffuse*(
                        (map_tmp[i+1][j][k]-2*map_tmp[i][j][k]+map_tmp[i-1][j][k])/(dx*dx)
                    +(map_tmp[i][j+1][k]-2*map_tmp[i][j][k]+map_tmp[i][j-1][k])/(dy*dy)
                    +(map_tmp[i][j][k+1]-2*map_tmp[i][j][k]+map_tmp[i][j][k-1])/(dz*dz)
                        );
        }
    }
}

for(int i=0;i<tNX;i++){
    for(int j=0;j<tNY;j++){
        for(int k=0;k<tNZ;k++){
            //alpha blending
            //weight‚Í0.0‚Æ1.0‚µ‚©Žg‚í‚È‚¢—\’è
            niche_value[i][j][k]=niche_gen_weight[i][j][k]*niche_gen[i][j][k]+(1.0-niche_gen_weight[i][j][k])*niche_value[i][j][k];

        }
    }
}
/*
*/
}
