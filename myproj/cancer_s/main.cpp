#include <iostream>
#include <memory>
#include "define.h"
#include "field.h"
//using namespace std;

int main(int argc, char *argv[])
{
    printf("oaioo");

    Field* fld=new Field();

    fld->define_niche_gen([](auto& gen,auto& weight){
for(int i=0;i<tNX;i++){
    for(int j=0;j<tNY;j++){
        for(int k=0;k<tNZ;k++){
            gen[i][j][k]=0;
            weight[i][j][k]=0;
        }
    }
}
for(int i=0;i<tNX;i++){
    for(int j=0;j<tNY;j++){
        gen[i][j][1]=1.0;
        weight[i][j][1]=1.0;
    }
}
    });

for(int i=0;i<100000;i++){
   std::cout << i << std::endl;
    fld->calc_niche_diffusion();
    std::cout << fld->niche_value[25][25][98] << std::endl;
}


    std::cout << "Hello World!" << std::endl;
    return 0;
}
