#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <algorithm>
static constexpr auto CELL_DATA_DIR = "output";
static constexpr auto IMG_OUT_DIR = "POV_out";
static constexpr auto TMP_FILE_NAME_PREF = "pov_tmp";
static constexpr double UMAX = 0.65;
static constexpr double UMIN = 0.1;

static constexpr double FMAX = 1.;
static constexpr double FMIN = 0.;


int MALIGNANT = 0;
int DISP_MODE = 0;
int DISP_MUSUME = 0;
int DISP_ALIVE = 0;
int DISP_DEAD = 0;
enum CELL_STATE:int {
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
std::string random_string(size_t length)
{
    static auto randchar = []() -> char
    {
        static const char charset[] =
            "0123456789"
            "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
            "abcdefghijklmnopqrstuvwxyz";
        const size_t max_index = (sizeof(charset) - 1);
        return charset[rand() % max_index];
    };
    std::string str(length, 0);
    std::generate_n(str.begin(), length, randchar);
    return str;
}
double ca_value_conv(double ca) {
    return ca < UMIN ? 0. :
        (ca > UMAX ? 1. :
            (ca - UMIN) / (UMAX - UMIN));
}

double fat_value_conv(double fat) {
    return fat < FMIN ? 0. :
        (fat > FMAX ? 1. :
            (fat - FMIN) / (FMAX - FMIN));
}

double infat_value_conv(double infat) {
    return infat < FMIN ? 0. :
        (infat > FMAX ? 1. :
            (infat - FMIN) / (FMAX - FMIN));
}

void setcolor_white(double *rgb)
{
    rgb[0] = 1.0;
    rgb[1] = 1.0;
    rgb[2] = 1.0;
}


void setcolor_reduction(double *rgb)
{
    rgb[0] *= 2.0;
    rgb[1] *= 2.0;
    rgb[2] *= 2.0;
}

void setcolor_musume(int div, double *rgb, int tb)
{
    double amp;
    if (tb < MALIGNANT) {
        rgb[0] = 1.0;
        rgb[1] = 0.7;
        rgb[2] = 0.7;
    }
    else {
        rgb[0] = 0.7;
        rgb[1] = 1.0;
        rgb[2] = 0.7;
    }

    if (div == 0) {
        rgb[0] += 0.2;
        rgb[1] += 0.2;
        rgb[2] += 0.2;
    }

}

void setcolor_fix(double *rgb,  int tb)
{
    if (tb < MALIGNANT) {
        rgb[0] = 0.8;
        rgb[1] = 0.0;
        rgb[2] = 0.1;
    }
    else {
        rgb[0] = 0.0;
        rgb[1] = 0.8;
        rgb[2] = 0.1;
    }
}

void setcolor_memb(double *rgb)
{
    rgb[0] = 0.5;
    rgb[1] = 0.5;
    rgb[2] = 0.5;
}
void setcolor_der(double *rgb)
{
    rgb[0] = 0.2;
    rgb[1] = 0.2;
    rgb[2] = 0.2;
}

void setcolor_gr(double *rgb, int tb)
{
    rgb[0] = 0.0;
    rgb[1] = 0.8;
    rgb[2] = 0.1;
}

void setcolor_ggr(double *rgb)
{
    rgb[0] = 0.;
    rgb[1] = 0.;
    rgb[2] = 1.;
}

void setcolor_ggr2(double *rgb)
{
    /*
    rgb[0] = 168. / 255.;
    rgb[1] = 87. / 255.;
    rgb[2] = 168. / 255. ;
    */
    rgb[0] = 0.3;
    rgb[1] = 0.3;
    rgb[2] = 0.3;

}

void setcolor_ggr3(double *rgb)
{

    rgb[0] = 247. / 255.;
    rgb[1] = 171. / 255.;
    rgb[2] = 166. / 255.;

}

void setcolor_pk(double *rgb)
{

    rgb[0] = 1.0;
    rgb[1] = 0.6;
    rgb[2] = 0.8;
    /*
    rgb[0] = 0.6;
    rgb[1] = 0.6;
    rgb[2] = 0.6;
    */
}
void setcolor_red(double *rgb)
{
    rgb[0] = 1.0;
    rgb[1] = 0.0;
    rgb[2] = 0.0;
}
void setcolor_pkwhite(double *rgb)
{

    rgb[0] = 1.0;
    rgb[1] = 0.8;
    rgb[2] = 1.0;
    /*
    rgb[0] = 0.6;
    rgb[1] = 0.6;
    rgb[2] = 0.6;
    */
}
void setcolor_M(double color, double *rgb, int tb)
{
    if (tb < MALIGNANT) {
        if (color < 0.0) {
            rgb[0] = 1.0;
            rgb[1] = 0.5;
            rgb[2] = 0.5;
        }
        else {
            rgb[0] = 1.0;
            rgb[1] = 0.0;
            rgb[2] = 0.0;
        }
    }
    else {
        if (color < 0.0) {
            rgb[0] = 1.0;
            rgb[1] = 0.6;
            rgb[2] = 0.8;
        }
        else if (color >= 1.0) {
            rgb[0] = 1.0;
            rgb[1] = 0.0;
            rgb[2] = 0.0;
        }
        else if ((color >= 0.0) && (color < 0.250)) {
            rgb[0] = 0.5;
            rgb[1] = color;
            rgb[2] = 0.5; //1.0
        }
        else if ((color >= 0.250) && (color < 0.50)) {
            rgb[0] = 0.0;
            rgb[1] = 1.0;
            rgb[2] = 2.0 - 4.0*color;
        }
        else if ((color >= 0.50) && (color < 0.750)) {
            rgb[0] = 4.0*color - 2.0;
            rgb[1] = 1.0;
            rgb[2] = 0.0;
        }
        else if ((color >= 0.750) && (color < 1.0)) {
            rgb[0] = 1.0;
            rgb[1] = 4.0 - 4.0*color;
            rgb[2] = 0.0;
        }
    }
}

void setcolor(double color, double *rgb)
{
    if (color < 0.0) {
        rgb[0] = 0.0;
        rgb[1] = 0.0;
        rgb[2] = 1.0;
    }
    else if (color >= 1.0) {
        rgb[0] = 1.0;
        rgb[1] = 0.0;
        rgb[2] = 0.0;
    }
    else if ((color >= 0.0) && (color < 0.250)) {
        rgb[0] = 0.0;
        rgb[1] = 4.0*color;
        rgb[2] = 1.0;
    }
    else if ((color >= 0.250) && (color < 0.50)) {
        rgb[0] = 0.0;
        rgb[1] = 1.0;
        rgb[2] = 2.0 - 4.0*color;
    }
    else if ((color >= 0.50) && (color < 0.750)) {
        rgb[0] = 4.0*color - 2.0;
        rgb[1] = 1.0;
        rgb[2] = 0.0;
    }
    else if ((color >= 0.750) && (color < 1.0)) {
        rgb[0] = 1.0;
        rgb[1] = 4.0 - 4.0*color;
        rgb[2] = 0.0;
    }
}
void setcolor2(double color, double *rgb)   //siro kuro
{
    if (color < 0.0) {
        rgb[0] = 1.0;
        rgb[1] = 1.0;
        rgb[2] = 1.0;
    }
    else if (color >= 1.0) {
        rgb[0] = 0.0;
        rgb[1] = 0.0;
        rgb[2] = 0.0;
    }
    else {
        rgb[0] = 1.0 - color;
        rgb[1] = 1.0 - color;
        rgb[2] = 1.0 - color;
    }
}
void setcolor3(double color, double *rgb)   //white blue
{
    if (color < 0.0) {
        rgb[0] = 1.0;
        rgb[1] = 1.0;
        rgb[2] = 1.0;
    }
    else if (color >= 2.0) {
        rgb[0] = 0.0;
        rgb[1] = 0.0;
        rgb[2] = 1.0;
    }
    else {
        rgb[0] = 1.0 - color;
        rgb[1] = 1.0 - color;
        rgb[2] = 1.0;
    }
}
void setcolor4(double color, double *rgb)   //white purple
{
    if (color < 0.0) {
        rgb[0] = 1.0;
        rgb[1] = 1.0;
        rgb[2] = 1.0;
    }
    else if (color >= 1.0) {
        rgb[0] = 0.6;
        rgb[1] = 0.2;
        rgb[2] = 0.8;
    }
    else {
        rgb[0] = 1.0 - 0.4*color;
        rgb[1] = 1.0 - 0.8*color;
        rgb[2] = 1.0 - 0.2*color;
    }
}

void calc_color(CELL_STATE state,double ca,double fat,double infat,int divtime,int tb,int touch,double *rgb) {
    if (state == MEMB) {
        if (DISP_MODE == 0) {
            setcolor_pkwhite(rgb);
        }
        else {
            setcolor_memb(rgb);
        }

        return;
    }

    if (state == FIX) {
        if (DISP_MODE == 0) {
            setcolor_M(0, rgb, tb);
        }
        else {
            setcolor_fix(rgb, tb);
        }
        return;
    }

    if (state == MUSUME) {
        if (DISP_MODE == 0) {
            setcolor_M(0, rgb, tb);
            return;
        }

        if (DISP_MUSUME == 1) {
            switch (DISP_MODE) {
            case 1:
                setcolor(ca_value_conv(ca), rgb);
                break;
            case 2:
                setcolor(fat_value_conv(fat), rgb);
                break;
            case 3:
                setcolor(infat_value_conv(infat), rgb);
                break;
            default:
                std::cout << "invalid param comb 1" << std::endl;
                exit(1);
                break;
            }
        }
        else if(DISP_MUSUME==2) {
            setcolor_musume(divtime, rgb, tb);
        }
        return;
    }

    if (state == ALIVE) {
        if (DISP_MODE == 0) {
            setcolor_M(0.122, rgb, tb);
            return;
        }

        if (DISP_ALIVE == 1 || (DISP_ALIVE==2&&touch==1)) {
            switch (DISP_MODE) {
            case 1:
                setcolor(ca_value_conv(ca), rgb);
                break;
            case 2:
                setcolor(fat_value_conv(fat), rgb);
                break;
            case 3:
                setcolor(infat_value_conv(infat), rgb);
                break;
            default:
                std::cout << "invalid param comb 2" << std::endl;
                exit(1);
                break;
            }
        }
        return;
    }
    

    if (state == DEAD) {
        switch (DISP_MODE) {
        case 0:
            setcolor_pk(rgb);
            break;
        case 1:
            setcolor_white(rgb);
            break;
        case 2:
            setcolor(fat_value_conv(fat), rgb);
            break;
        case 3:
            setcolor(infat_value_conv(infat), rgb);
            break;
        default:
            std::cout << "invalid param comb 3" << std::endl;
            exit(1);
            break;
        }
        return;
    }
}
bool should_render(CELL_STATE state){
	if(state==MEMB)return true;
	if(state==FIX)return true;
	if(state==MUSUME&&DISP_MUSUME!=0)return true;
	if(state==ALIVE&&DISP_ALIVE!=0)return true;
	if(state==DEAD&&DISP_DEAD!=0)return true;
	return false;
}
void set_pov_template(std::ofstream& ofs) {
    ofs <<
        "#version 3.6;\n"
        "#include \"colors.inc\"\n"
        "#include \"metals.inc\"\n"
        "#include \"textures.inc\"\n"
        "global_settings{max_trace_level 16 }\n"
        "camera{\n"
        "location <50,-300,80>\n"
        "sky z\n"
        "right 0.24*x*image_width / image_height\n"
        "up 0.24*z\n"
        "look_at <50,25,40>}\n"
        "background{ rgb 0 }\n"
      //"light_source{ <-8,-20,50> color rgb <0.2,0.2,0.2> shadowless }\n"
      // "light_source{ <108,-20,50> color rgb <0.2,0.2,0.2> shadowless }\n"
        "light_source{ <50,-100,50> color rgb <0.45,0.45,0.45>  shadowless}\n"
      //"light_source{ <50,25,-200> color rgb <0.1,0.1,0.1> shadowless}\n"
      //  "light_source{ <50,25,200> color rgb <0.1,0.1,0.1> shadowless}\n"
;
}
int main(int argc, char* argv[]) {
    int DISP_ORIGIN_BEGIN=0;
    int DISP_ORIGIN_END=0;

    if (argc < 8) {
        exit(1);
    }
    printf("%d\n", argc);
    if (argc>9) {
        if (argc == 10) {
            printf("NO ORIGIN_END SPECIFIED\n"); exit(1);
        }
        DISP_ORIGIN_BEGIN = atoi(argv[9]);
        DISP_ORIGIN_END = atoi(argv[10]);
    }


    std::string datadir = argv[1];
    MALIGNANT = atoi(argv[2]);

    int t0 = atoi(argv[3]);
    int tmax = atoi(argv[4]);

    DISP_MODE = atoi(argv[5]);
    DISP_MUSUME = atoi(argv[6]);
    DISP_ALIVE = atoi(argv[7]);
    DISP_DEAD = atoi(argv[8]);
    bool DISP_FLAGS[10] = { true };
    if (DISP_MUSUME == 0)DISP_FLAGS[MUSUME] = false;
    if (DISP_ALIVE == 0)DISP_FLAGS[ALIVE] = false;
    if (DISP_DEAD == 0)DISP_FLAGS[DEAD] = false;
    //printf("x0=%f y0=%f DISP_MODE=%d\n", xo, yo, DISP_MODE);
system(("mkdir "+std::string(IMG_OUT_DIR)).c_str());
    for (int i = t0; i < tmax; i++) {
        std::string datapath = std::string(CELL_DATA_DIR) + "/" + std::to_string(i);
        std::ifstream cdata(datapath);
        if (!cdata) {
            std::cout << "An error has occured when opening the file:" << datapath << std::endl;
            exit(1);
        }
        std::string tmp_file = std::string(TMP_FILE_NAME_PREF) + random_string(8) + "~";
        std::ofstream povdata(tmp_file,std::ios_base::out|std::ios_base::trunc);
        if (!povdata) {
            std::cout << "An error has occured when creating a temporary file:" << tmp_file << std::endl;
            exit(1);
        }
        set_pov_template(povdata);

       // povdata << "union{" << std::endl;
        std::string line;
        while (std::getline(cdata, line)) {
            CELL_STATE state = UNUSED;
            double rad = 0.0;
            double agek = 0.0;
            double x = 0.0, y = 0.0, z = 0.0;
            double ca = 0.0;
            int divtime = -1;
            double fat = 0.0;
            double infat = 0;
            int touch = 0;
            int tb = -1;
            sscanf(line.c_str(), "%*d %d %lf %*lf %lf %*lf %lf %lf %lf %lf %d %lf %lf %d %*lf %*d %d",
                (int*)&state, &rad, &agek, &x, &y, &z, &ca, &divtime, &fat, &infat, &touch, &tb);
		if(!should_render(state))continue;
            double rgb[3];
            calc_color(state, ca, fat, infat, divtime, tb, touch, rgb);
            povdata << "sphere{<" 
                << x << "," << y << "," << z << ">," 
                << rad << " "
                <<"texture{pigment{color rgb<"<<rgb[0]<<","<<rgb[1]<<","<<rgb[2]<<">} finish{diffuse 2.5}}}" << std::endl;

            
        }
        //povdata << "}" << std::endl;
std::cout<<"Rendering "<<i<<" ..."<<std::endl;
        system(("povray +I"+tmp_file+" +H400 +W800 +O"+IMG_OUT_DIR+"/"+std::to_string(i)+" >/dev/null 2>&1").c_str());
        system(("rm -f " + tmp_file).c_str());
std::cout<<"done."<<std::endl;
    }
}
