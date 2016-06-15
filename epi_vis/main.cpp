#include <iostream>
#include <fstream>
#include <string>
#include <voro/voro++.hh>
using namespace std;

static constexpr double LX=50.0;
static constexpr double LY=50.0;
static constexpr double LZ=100.0;

static constexpr int comp_nx=8,comp_ny=8,comp_nz=8;

enum CELL_STATE {
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
int main(int argc, char *argv[])
{
    if(argc<=1){
        cout << "No cell data specified." << endl;
        exit(1);
    }
    ifstream ifs(argv[1]);
    if(!ifs){
        cout << "Opening cell data failed. Filename:"<<argv[1] << endl;
        exit(1);
    }

    voro::container con(0,LX,0,LY,0,LZ,comp_nx,comp_ny,comp_nz,true,true,false,8);
    string str;
    int count=0;
    while(getline(ifs,str)){
        int index=0;
        int state=0;
        double x,y,z;
        sscanf(str.data(),"%d %d %*f %*f %*f %*f %lf %lf %lf %*f %*d %*f %*f %*d %*f %*d %*d",&index,&state,&x,&y,&z);
        if(state==MUSUME||state==ALIVE||state==DEAD||state==FIX){
        con.put(count++,x,y,z);
        }
    }
    con.draw_cells_pov("cell_voro_test.pov");
    //cout << "Hello World!" << endl;
    return 0;
}
