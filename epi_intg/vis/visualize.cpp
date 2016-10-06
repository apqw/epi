#include "../global.h"
#include "visualize.h"
#include "VisParams.h"

#include "../calc/cpu2/CellManager.h"
//#include "../opencsg/opencsg.h"
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
//#include <gl/osmesa.h>
#include "../lodepng/lodepng.h"
#include <algorithm>
#include "../fsys.h"
#include "vis_global.h"
#include "color.h"

static constexpr int SPHERE_DIV = 16;
static bool initialized = false;
static VisParams __dummy_vp=VisParams();
static const GLdouble fovy = 30.0;
static const GLdouble znear = 1.0;
static GLdouble zfar=-1.0; //uninitialized
static GLdouble distance = -1.0;
static const GLdouble twist = 0.0;
static std::vector<uint8_t> raw_buf;
static double OFFSETX;
static double OFFSETY;
static std::string output_path="";
VisParams& vp = __dummy_vp;


static void set_initial_var(const VisParams&);

static void set_view();
static void set_polarview();

static void initLighting();
static void set_initial_view();
static void init_internal();
static void fix_coord(double inx, double iny, double inz, double*o1, double*o2, double*o3);
static void flipYPixel(std::vector<uint8_t>& pix, int w, int h, int ch);

static void clear_gl_buffer();
static void begin_draw_setting();

static void end_draw_setting();

static void output_img_from_buffer(int index);







void init_visualizer(int* argcp, char** argv, const VisParams& vp) {
    // glutInit(argcp, argv);
    glutInit(argcp, argv);
    set_initial_var(vp);
    initialized = true;
}


static void set_view() {
    GLdouble aspect = (double)vp.width / (double)vp.height;


    glViewport(0, 0, vp.width, vp.height);    /* hyouzi hanni setttei */
    glMatrixMode(GL_PROJECTION);        /* touei houhou settei */
    glLoadIdentity();                   /* zahyouhennkanngyouretu syokika */
    gluPerspective(fovy, aspect, znear, zfar);
    glMatrixMode(GL_MODELVIEW);
}

static void set_polarview() {
    glTranslatef(0.0, 0.0, -distance);
    glRotatef(-twist, 0.0, 0.0, 1.0);
    glRotatef(-vp.elevation, 1.0, 0.0, 0.0);
    glRotatef(-vp.azimuth, 0.0, 1.0, 0.0);
}

static void set_initial_var(const VisParams&_vp) {
    vp = _vp;
    zfar = 4.0*pm->LX;
    distance = 0.6*zfar;
    raw_buf.resize(vp.width*vp.height * 4);
    OFFSETX = pm->LX / 2.0;
    OFFSETY = pm->LY / 2.0;
    output_path = vp.output + "/";
}



static void initLighting()
{
    static const     float  diffuse[] = { 1.0f, 1.0f, 1.0f, 1.0f };
    static const     float  specular[] = { 0.8f, 0.8f, 0.8f, 1.0f };
static const float  ambient[] = { 0.1f, 0.1f, 0.1f, 1.0f };

    glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);     /*gl_dif(kakusannkou): */
    glLightfv(GL_LIGHT0, GL_SPECULAR, specular);   /*gl_sp(kyoumennkou):(1,1,1,1)*/
    glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);    /*gl_am(kannkyoukou):(0,0,0,1)*/
    glEnable(GL_LIGHT0);

    static float light[] = { 1.0f, 1.0f, 5.0f, 0.0f };
    glLightfv(GL_LIGHT0, GL_POSITION, light);
}


static void set_initial_view() {
    set_view();
    set_polarview();
}


static void init_internal() {
    glutInitWindowPosition(0, 0);
    glutInitWindowSize(vp.width, vp.height);
    glutInitDisplayMode(/* GLUT_SINGLE */ GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
    glutCreateWindow("epid_vis");
    glClearColor(0.0, 0.0, 0.0, 1.0);
    glShadeModel(GL_SMOOTH);
    initLighting();
    glReadBuffer(GL_FRONT);

    

    make_dir(output_path.c_str());


    set_initial_view();
}



static void fix_coord(double inx, double iny, double inz, double*o1, double*o2, double*o3) {
    *o1 = inx - vp.X0-OFFSETX;
    *o2 = -inz + pm->LZ / 4.0 - vp.z_level;
    *o3 = iny - vp.Y0-OFFSETY;

}
static void flipYPixel(std::vector<uint8_t>& pix,int w,int h,int ch) {
    for (int i = 0; i < h / 2; i++) {
        std::swap_ranges(pix.begin() + i*w*ch, pix.begin() + (i+1)*w*ch, pix.begin() + (h-i-1)*w*ch);
    }
}

static bool should_draw(const CellTempStruct&cts) {

    switch (cts.state) {
    case BLANK:case DER:
        return false;
        break;
    case MEMB:
        return true;
        break;
    case FIX:
        return true;
        break;
    case MUSUME:
        return vp.disp_musume != 0;
        break;
    case ALIVE:
        if (vp.disp_alive == 2) {
            return cts.touch == 1;
        }
        else {
            return vp.disp_alive != 0;
        }
        break;
    case DEAD:
        return vp.disp_dead != 0;
        break;
    default:
        return false;
        break;
    }
}


struct OnCellLoadVis {
    //const VisParams& vp;
    static const float diffuse[4];// = { 1.0, 1.0, 1.0, 1.0 };
    static const float specular[4];// = { 0.3, 0.3, 0.3, 1.0 };
    static const float ambient[4];// = { 0.3, 0.3, 0.3, 1.0 };
    OnCellLoadVis() {
    }
    void operator()(CellManager&cman, const CellTempStruct&cts) {

        if (!should_draw(cts))return;
        

        std::vector<float> color = get_color_rgb(cts);
        color.push_back(1.0f); //if you need Alpha
        glMaterialfv(GL_FRONT, GL_DIFFUSE, &color[0]);
        glMaterialfv(GL_FRONT, GL_SPECULAR, specular);
        glMaterialfv(GL_FRONT, GL_AMBIENT, ambient);
        glMaterialf(GL_FRONT, GL_SHININESS, 128.0);

        double ox, oy, oz;
        fix_coord( cts.x, cts.y, cts.z, &ox, &oy, &oz);
        glTranslatef(ox, oy, oz); //xzy

        glutSolidSphere(cts.rad, SPHERE_DIV, SPHERE_DIV);

        glTranslatef(-ox, -oy, -oz);

        //cman.test_realloc();
    }
};

const float OnCellLoadVis::diffuse[4] = { 1.0f, 1.0f, 1.0f, 1.0f };
const float OnCellLoadVis::specular[4] = { 0.3f, 0.3f, 0.3f, 1.0f };
const float OnCellLoadVis::ambient[4] = { 0.3f, 0.3f, 0.3f, 1.0f };
static void clear_gl_buffer() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0.0, 0.0, 0.0, 1.0);
}

static void begin_draw_setting() {
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
}

static void end_draw_setting() {
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);
}

static void output_img_from_buffer(int index) {
    glReadPixels(0, 0, vp.width, vp.height, GL_RGBA, GL_UNSIGNED_BYTE, &raw_buf[0]);
    flipYPixel(raw_buf, vp.width, vp.height, 4);
    lodepng::encode(output_path + std::to_string(index) + ".png", raw_buf, vp.width, vp.height);
}


void visualize() {
    if (!initialized) {
        throw std::logic_error("Visualizer not initialized.");
    }

    
    init_internal();

    CellManager cman;
    
    auto on = OnCellLoadVis();
    std::cout<<"Visualize:\nGenerating images..."<<std::endl;
    for (int i = vp.start; i <= vp.end; i++) {
        clear_gl_buffer();

        begin_draw_setting();
        
        glPushMatrix();
        cman.load(vp.datadir+"/"+std::to_string(i), on);
        
        glPopMatrix();
        end_draw_setting();

        glutSwapBuffers();

        output_img_from_buffer(i);
        std::cout<<"Done "<<i<<std::endl;
    }
}
