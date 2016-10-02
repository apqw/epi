#pragma once
#include "../Params.h"
#include "../define.h"
class VisParams:public Params
{
public:
    std::string datadir;
    std::string output;
    unsigned int MALIGNANT;
    
    int start;
    int end;
    real X0, Y0;
    real elevation;
    real azimuth;
    real z_level;
    int mode;
    int disp_musume;
    int disp_alive;
    int disp_dead;
    unsigned int width, height;

    enum DISP_MODE:int {
        MEDICAL=0,CA2P=1,EX_FAT=2,IN_FAT=3
    };

    VisParams();
    VisParams(const std::string& paramfile);
    void init();

    //void load(const std::string& paramfile);
};

