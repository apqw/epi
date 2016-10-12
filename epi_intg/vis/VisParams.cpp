#include "VisParams.h"



VisParams::VisParams()
{
    s_Ctor();
}


VisParams::VisParams(const std::string& path)
{
    s_Ctor(path);
}

void VisParams::init() {
    piset = {
        gp1(datadir),
        gp1(output),
        gp1(MALIGNANT),
        gp1(start),
        gp1(end),
        gp1(X0),
        gp1(Y0),
        gp1(elevation),
        gp1(azimuth),
        gp1(z_level),
        gp1(mode),
        gp1(disp_musume),
        gp1(disp_alive),
        gp1(disp_dead),
        gp1(width),
        gp1(height),
gp1(color_height_min),
gp1(color_height_max),
gp1(TEMP_DIRTY1)
    };
    datadir = "";
    output = "3Dimage_tmp";
    MALIGNANT = 0;
    start = 0;
    end = 0;
    X0 = 0;
    Y0 = 0;
    elevation = 0.0;
    azimuth = 175;
    z_level = 0;
    mode = 0;
    disp_musume = 1;
    disp_alive = 1;
    disp_dead = 1;
    width = 320;
    height = 320;
color_height_min=0.0;
color_height_max=50;
TEMP_DIRTY1=false;
}
