/*
  ¥Õ¥¡¥¤¥ë¤«¤éÇÛÃÖ¤ÈÃÍ¤òÆÉ¤ß¹þ¤ß¡¤É½¼¨¤¹¤ë¡¥
  »þ´ÖÊÑ²½¤Ë¤âÂÐ±þ
  ¥º¡¼¥à¤ÎÊÑ¹¹
  JPEG¥Õ¥¡¥¤¥ë¤Ë½ÐÎÏ¤¹¤ë
  È¾·Â¤ò¥Õ¥¡¥¤¥ë¤«¤éÆÉ¤ß¹þ¤à
  ³«»ÏÈÖ¹æ¡¤½ªÎ»ÈÖ¹æ¤ò¥³¥Þ¥ó¥É¥é¥¤¥ó¤«¤éÆþÎÏ²Ä
  (¥Õ¥¡¥¤¥ë¤«¤éÆÉ¤ß¹þ¤ß¤Ê¤¬¤éÉÁ²è¤¹¤ë)
  zÊý¸þ¤ò¸µ¤ËÌá¤¹
  ------------- OpenGL1-R-3.c ¤È OpenGL2-R-3.c ¤òÅý¹ç -----------------
*/

/*  c7-1.c   Copyright (c) 2003 by T. HAYASHI and K. KATO  */
/*                                    All rights reserved  */
#include <stdio.h>
#include <stdlib.h>
#include <GL/glut.h>
#include <GL/glu.h>
#include <GL/gl.h>
#include <math.h>
#include <GL/osmesa.h>
#include <unistd.h>
//#include "param.h"

#define SAVE 1

#define KEY_ESC 27

/* state list */
#define ALIVE 0
#define DEAD  1
#define DISA  2
#define ATP   3
#define FIX   4
#define BLANK 5
#define DER 6
#define MUSUME 7
#define MEMB 9

#define NN 30000
#define N NN

#define MAXTIME 5000
#define sample 10 /* µå¤Î¤Ê¤á¤é¤«¤µ */
#define UMAX 0.65 
#define UMIN 0.1
#define FMAX 1.0
#define FMIN 0.0
#define Ne 200

/*
define DISTANCE 200.0
#define ELEVATION (-90.)
#define AZIMUTH 0.
*/
double ZLEVEL;
int time = 0;
int Width =  400;
int Height = 400;
double view_near = 100.;
double fovy = 30.; //15.
double znear = 1.;

double LX=50;
double LY=50;
double LZ=100; // param.h
int key;

int MALIGNANT;
int DISP_MODE;
int DISP_MUSUME;
int DISP_ALIVE;
int DISP_DEAD;

double ELEVATION; //170.0//90.0//(90.0)

double AZIMUTH; //0.0//90.0

#define GLX 0
#define GLY 750
/* #include"para_GL.h" */
#define OUTPUT_IMAGE "3Dimage_tmp" //画像の保存先

void polarview( void );
void resetview( void );
void myReshape(int, int);
void save_file(FILE *f, void *buffer);
void change_image_format(char *in, char *out);
void setcolor_white(double rgb[3]);
void setcolor_musume(int div, double rgb[3], int tb);

void setcolor_gr(double rgb[3], int tb);
void setcolor_der(double rgb[3]);
void setcolor_memb(double *rgb);
void setcolor_fix(double rgb[3], int tb);
void setcolor_ggr(double rgb[3]);
void setcolor_ggr2(double rgb[3]);
void setcolor_ggr3(double rgb[3]);
void Drawsycle(double s,double rad);

void setcolor_pk(double rgb[3]);
void setcolor_red(double rgb[3]);
void setcolor_pkwhite(double *rgb);
void setcolor(double color, double rgb[3]);
void setcolor_M(double color, double rgb[3], int tb);
void setcolor2(double color, double rgb[3]);
void setcolor3(double color, double rgb[3]);
void setcolor4(double color, double rgb[3]);
void setorg(void);
double frand(void);
double dist(double x[3], double y[3]);
double inner(double x[3], double y[3]);
void input(void);

float distance, twist, elevation, azimuth;
int xBegin, yBegin;
int mButton;
float theta =  15.;
unsigned char revolveFlag = GL_FALSE;

int t0, tmax;
int state[NN];
double xyz[NN][3], value[NN], agek[NN], fat[NN], ca[NN];
int tb[NN];
double rad[NN];
double org[3];
double org_e[3];
char *datadir;

double xo, yo;

int divtime[NN], touch[NN];

char done = 0;

GLubyte *buffer;

double zfar;
double DISTANCE;

void display(void)
{
  int i,j;
  int n = NN;
  double x, y, z;
  double r, g, b;
  float diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
  float specular[] = { 0.3, 0.3, 0.3, 1.0 };
  float ambient[] = { 0.3, 0.3, 0.3, 1.0 };
  float light[] = { 1.0, 1.0, 5.0, 0.0 };
  double tmpcolor[3];
  char ppmfile[256], jpegfile[256];
  FILE *fp;
  
  glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
  glClearColor (0.0, 0.0, 0.0, 1.0); // background color
  //  glClearColor (1.0, 1.0, 1.0, 1.0); // background color
  glEnable( GL_DEPTH_TEST );
  glEnable( GL_LIGHTING );
  
  glPushMatrix();           /* gennzai zahyoukei hozonn*/
  
  glLightfv( GL_LIGHT0, GL_POSITION, light );   /*GL_POSITION(syoumeiiti):(0,0,1,0)*/

  myReshape(Width, Height);

  polarview();
  
  for(i = 0; i < n; i++) {
    if (state[i] == BLANK) break;

    else if (state[i] == MEMB) {
      if (DISP_MODE==0)
	setcolor_pkwhite(tmpcolor);
      else
	setcolor_memb(tmpcolor);
    }
    else if(state[i] == FIX){
      if (DISP_MODE==0) // medical
	setcolor_M(-1.0, tmpcolor, tb[i]);
      else
	setcolor_fix(tmpcolor, tb[i]);
    }
    else if(state[i] == MUSUME && DISP_MUSUME != 0) {
      if (DISP_MODE==0) 
 	setcolor_M(0.0, tmpcolor, tb[i]);
      else if (DISP_MUSUME==1) {
	if (DISP_MODE==1)
	  setcolor(ca[i] < UMIN ? 0. :
		   (ca[i] > UMAX ? 1. :
		    (ca[i] - UMIN) / (UMAX - UMIN)), tmpcolor);
	else if (DISP_MODE==2)
	  setcolor(fat[i] < FMIN ? 0. :
		   (fat[i] > FMAX ? 1. :
		    (fat[i] - FMIN) / (FMAX - FMIN)), tmpcolor);
      }
      else if (DISP_MUSUME==2)
	setcolor_musume(divtime[i], tmpcolor, tb[i]);
    }
    else if (state[i] == ALIVE && DISP_ALIVE != 0) {
      if (DISP_MODE==0) 
 	setcolor_M(0.122, tmpcolor, tb[i]);
      else if (DISP_ALIVE==1 || (DISP_ALIVE==2 && touch[i] == 1)) {
	if (DISP_MODE==1)
	  setcolor(ca[i] < UMIN ? 0. :
		   (ca[i] > UMAX ? 1. :
		    (ca[i] - UMIN) / (UMAX - UMIN)), tmpcolor);
	else if (DISP_MODE==2)
	  setcolor(fat[i] < FMIN ? 0. :
		   (fat[i] > FMAX ? 1. :
		    (fat[i] - FMIN) / (FMAX - FMIN)), tmpcolor);
	
      }
      else continue;
    }
    else if (state[i] == DEAD && DISP_DEAD) {
      if (DISP_MODE==0)
	setcolor_pk(tmpcolor);
      else if (DISP_MODE==1)
	setcolor_white(tmpcolor);
      else if (DISP_MODE==2)
	  setcolor(fat[i] < FMIN ? 0. :
		   (fat[i] > FMAX ? 1. :
		    (fat[i] - FMIN) / (FMAX - FMIN)), tmpcolor);
    }
    else continue;

    diffuse[0] = tmpcolor[0];
    diffuse[1] = tmpcolor[1];
    diffuse[2] = tmpcolor[2];
    glMaterialfv( GL_FRONT, GL_DIFFUSE, diffuse );      /*zennmenn kakusannkou*/
    glMaterialfv( GL_FRONT, GL_SPECULAR, specular );    /*kyoumennkolu*/
    glMaterialfv( GL_FRONT, GL_AMBIENT, ambient );      /*kannkyoukou */
    glMaterialf( GL_FRONT, GL_SHININESS, 128.0 );        /*shin:kyoumennkou sisuu*/
    
    glTranslatef(xyz[i][0], xyz[i][1], xyz[i][2]);
    
    glutSolidSphere(rad[i], sample, sample);  
    
    glTranslatef(-xyz[i][0], -xyz[i][1], -xyz[i][2]);
  }
  glDisable( GL_LIGHTING );

  glDisable( GL_DEPTH_TEST );

  glPopMatrix();
  
  //glFlush();
  glutSwapBuffers();
  
  if(SAVE == 1) {
    buffer = malloc(Width * Height * 4);
    glReadPixels(0, 0, Width, Height, GL_RGBA, GL_UNSIGNED_BYTE, buffer);

    sprintf(ppmfile, "%s/c%05d.ppm",OUTPUT_IMAGE, time);
    sprintf(jpegfile, "%s/c%05d.jpg",OUTPUT_IMAGE, time);

    if((fp = fopen(ppmfile, "w")) == NULL) {
      printf("output file error\n");
      exit(1);
    }
    save_file(fp, buffer);
    change_image_format(ppmfile,jpegfile);
    remove(ppmfile);
    free(buffer);
    if(time == tmax) {
      exit(0);
    }
  } else {
    if(time == tmax && done == 0) {
      printf("End\n");
      done = 1;
    }
  }
  
  if(done == 0) {
    time += 1;
    input();
  }
  }
  
  void change_image_format(char *in, char *out)
{
  char string[300];
  sprintf(string, "convert %s %s", in, out);
  system(string);
}

void initLighting( void )
{
  float  diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
  float  specular[] = { 0.8, 0.8, 0.8, 1.0 };
  float  ambient[] = { 0.1, 0.1, 0.1, 1.0 };
  
  glLightfv( GL_LIGHT0, GL_DIFFUSE, diffuse );     /*gl_dif(kakusannkou): */
  glLightfv( GL_LIGHT0, GL_SPECULAR, specular );   /*gl_sp(kyoumennkou):(1,1,1,1)*/
  glLightfv( GL_LIGHT0, GL_AMBIENT, ambient );    /*gl_am(kannkyoukou):(0,0,0,1)*/
  glEnable( GL_LIGHT0 );
}

void input(void)
{
  int n = N, i;
  double x, y, z, ub, dumx, dumy;
  char str[256];
  FILE *fi;

  //sprintf(str, "%s/data%05d/%d", datadir, time, time);
  sprintf(str, "%s/%d", datadir, time);
  if((fi = fopen(str, "r")) == NULL) {
    printf("input file error\n");
    exit(1);
  }
  for(i = 0; i < n; i++) {
    fscanf(fi, "%*d %d %lf %*lf %lf %*lf %lf %lf %lf %lf %d %lf %*lf %d %*lf %*d %d", 
	   &state[i], &rad[i], &agek[i], &x, &y, &z, &ca[i], &divtime[i], &fat[i], &touch[i], &tb[i]);

    dumx = x - xo + LX/2.0;
    dumy = y - yo + LY/2.0;

    while (dumx > LX) {
      dumx -= LX;
    }
    while (dumx < 0.0) {
      dumx += LX;
    }
    while (dumy > LY) {
      dumy -= LY;
    }
    while (dumy < 0.0) {
      dumy += LY;
    }

    xyz[i][0] = dumx - LX / 2.;
    xyz[i][1] = - z + LZ/4 -ZLEVEL;
    xyz[i][2] = dumy - LY / 2.;
  }
  fclose(fi);
}

void idle(void)
{
  glutPostRedisplay();
}

void myKbd( unsigned char key, int x, int y )
{
  switch(key) {
  case KEY_ESC:
    exit( 0 );
  case 'q':
    exit(0);
  }
}

void myMouse( int button, int state, int x, int y )
{
  if (state == GLUT_DOWN) {
    switch(button) {
    case GLUT_LEFT_BUTTON:
      mButton = button;
      break;
    case GLUT_MIDDLE_BUTTON:
      revolveFlag = !revolveFlag;
      if( revolveFlag == GL_TRUE )
	glutIdleFunc( idle );
      else
	glutIdleFunc(NULL);
      break;
    case GLUT_RIGHT_BUTTON:
      mButton = button;
      break;
    }
    xBegin = x;
    yBegin = y;
  }
}


void myMotion( int x, int y )
{
  int xDisp, yDisp;
	
  xDisp = x - xBegin;
  yDisp = y - yBegin;

  switch (mButton) {
  case GLUT_LEFT_BUTTON:
    azimuth -= (float) xDisp/2.0;
    elevation -= (float) yDisp/2.0;
    break;
  case GLUT_RIGHT_BUTTON:
    //distance += (float) yDisp/4.;
    view_near -= (float) yDisp/4.;/* sisenn tikadukeru */
    break;
  }
  xBegin = x;
  yBegin = y;
  setorg();
  glutPostRedisplay();
}


void myInit (char *progname)
{
  glutInitWindowPosition(GLX, GLY);
  glutInitWindowSize( Width, Height );
  glutInitDisplayMode(/* GLUT_SINGLE */ GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
  glutCreateWindow(progname);
  glClearColor (0.0, 0.0, 0.0, 1.0);
  glutKeyboardFunc( myKbd );      /*gamen syoukyo settei */
  glutMouseFunc(myMouse);
  glutMotionFunc(myMotion);
  resetview();
  glShadeModel( GL_SMOOTH );
  initLighting();
  setorg();
  glReadBuffer(GL_FRONT);
}


void myReshape(int width, int height)
{
  float aspect = (float) width / (float) height;
  Width = width;
  Height = height;
  glViewport(0, 0, width, height);    /* hyouzi hanni setttei */
  glMatrixMode(GL_PROJECTION);        /* touei houhou settei */
  glLoadIdentity();                   /* zahyouhennkanngyouretu syokika */
  gluPerspective(fovy, aspect, znear, zfar);     
  glMatrixMode(GL_MODELVIEW);
}


void polarview( void )
{
  glTranslatef( 0.0, 0.0, -distance);
  glRotatef( -twist, 0.0, 0.0, 1.0);
  glRotatef( -elevation, 1.0, 0.0, 0.0);
  glRotatef( -azimuth, 0.0, 1.0, 0.0);
}


void resetview( void )
{
  distance = DISTANCE;
  twist = 0.;
  elevation = ELEVATION;
  azimuth = AZIMUTH;
}


int main(int argc, char** argv)
{
  int i, j;
  int n = N;
  double x, y, z;
  FILE *fi;
  char str[256];

  if(argc < 10) {
    exit(1);
  }

  datadir = argv[1];
  MALIGNANT = atoi(argv[2]);
  zfar = (4*LX);
  DISTANCE = 0.6*zfar;

  t0 = atoi(argv[3]);
  tmax = atoi(argv[4]);

  xo = atof(argv[5]);
  yo = atof(argv[6]);

  ELEVATION = atof(argv[7]);
  AZIMUTH = atof(argv[8]);
  ZLEVEL = atof(argv[9]);

  DISP_MODE = atoi(argv[10]);
  DISP_MUSUME = atoi(argv[11]);
  DISP_ALIVE = atoi(argv[12]);
  DISP_DEAD = atoi(argv[13]);

  printf("x0=%f y0=%f DISP_MODE=%d\n", xo, yo, DISP_MODE);

  sprintf(str, "mkdir %s", OUTPUT_IMAGE);
  system(str);

  time = t0;

  input();

  glutInit(&argc, argv);
  myInit(argv[0]);               /* syokisettei? */


  glutReshapeFunc(myReshape);    /*saizuhennkou kannsuu touroku*/

  glutDisplayFunc(display);      /*gamennkousinnyou kannsuu touroku*/
  glutIdleFunc(idle);
  glutMainLoop(); 
  return( 0 );
}

double frand(void)
{
  return (double)rand() / (double)RAND_MAX;
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

void setcolor_fix(double *rgb, int tb)
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
  rgb[2] = 166. / 255. ;
  
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
    if(color < 0.0){ 
      rgb[0] = 1.0;
      rgb[1] = 0.6;
      rgb[2] = 0.8;
    }else if(color >= 1.0){
      rgb[0] = 1.0;
      rgb[1] = 0.0;
      rgb[2] = 0.0;
    }else if( (color >= 0.0) && (color < 0.250)){
      rgb[0] = 0.5;
      rgb[1] = color;
      rgb[2] = 0.5; //1.0
    }else if( (color >= 0.250) && (color < 0.50)){
      rgb[0] = 0.0;
      rgb[1] = 1.0;
      rgb[2] = 2.0 - 4.0*color;
    }else if( (color >= 0.50) && (color < 0.750)){
      rgb[0] = 4.0*color - 2.0;
      rgb[1] = 1.0;
      rgb[2] = 0.0;
    }else if( (color >= 0.750) && (color < 1.0)){
      rgb[0] = 1.0;
      rgb[1] = 4.0 - 4.0*color;
      rgb[2] = 0.0;
    } 
  }
}

void setcolor(double color, double *rgb)
{
  if(color < 0.0){ 
    rgb[0] = 0.0;
    rgb[1] = 0.0;
    rgb[2] = 1.0;
  }else if(color >= 1.0){
    rgb[0] = 1.0;
    rgb[1] = 0.0;
    rgb[2] = 0.0;
  }else if( (color >= 0.0) && (color < 0.250)){
    rgb[0] = 0.0;
    rgb[1] = 4.0*color;
    rgb[2] = 1.0;
  }else if( (color >= 0.250) && (color < 0.50)){
    rgb[0] = 0.0;
    rgb[1] = 1.0;
    rgb[2] = 2.0 - 4.0*color;
  }else if( (color >= 0.50) && (color < 0.750)){
    rgb[0] = 4.0*color - 2.0;
    rgb[1] = 1.0;
    rgb[2] = 0.0;
  }else if( (color >= 0.750) && (color < 1.0)){
    rgb[0] = 1.0;
    rgb[1] = 4.0 - 4.0*color;
    rgb[2] = 0.0;
  } 
}
void setcolor2(double color, double *rgb)   //siro kuro
{
  if(color < 0.0){ 
    rgb[0] = 1.0;
    rgb[1] = 1.0;
    rgb[2] = 1.0;
  }else if(color >= 1.0){
    rgb[0] = 0.0;
    rgb[1] = 0.0;
    rgb[2] = 0.0;
  }else{
    rgb[0] =1.0-color;
    rgb[1] =1.0-color;
    rgb[2] =1.0-color;
  } 
}
void setcolor3(double color, double *rgb)   //white blue
{
  if(color < 0.0){ 
    rgb[0] = 1.0;
    rgb[1] = 1.0;
    rgb[2] = 1.0;
  }else if(color >= 2.0){
    rgb[0] = 0.0;
    rgb[1] = 0.0;
    rgb[2] = 1.0;
  }else{
    rgb[0] =1.0-color;	
    rgb[1] =1.0-color;
    rgb[2] =1.0;
  } 
}
void setcolor4(double color, double *rgb)   //white purple
{
  if(color < 0.0){ 
    rgb[0] = 1.0;
    rgb[1] = 1.0;
    rgb[2] = 1.0;
  }else if(color >= 1.0){
    rgb[0] = 0.6;
    rgb[1] = 0.2;
    rgb[2] = 0.8;
  }else{
    rgb[0] =1.0-0.4*color;	
    rgb[1] =1.0-0.8*color;
    rgb[2] =1.0-0.2*color;
  } 
}

double dist(double x[3], double y[3])
{
  return sqrt((x[0] - y[0]) * (x[0] - y[0])
	      + (x[1] - y[1]) * (x[1] - y[1])
	      + (x[2] - y[2]) * (x[2] - y[2]));
}

void setorg(void)
{
  double l;
  /*
  org[0] = view_near * cos(-elevation * M_PI / 180.) * sin(azimuth * M_PI / 180.);
  org[1] = view_near * sin(-elevation * M_PI / 180.);
  org[2] = view_near * cos(-elevation * M_PI / 180.) * cos(azimuth * M_PI / 180.);
  */
  org[0] = view_near * sin(azimuth * M_PI / 180.);
  org[1] = 0.;
  org[2] = view_near * cos(azimuth * M_PI / 180.);
  l = sqrt(org[0] * org[0] + org[1] * org[1] + org[2] * org[2]);
  org_e[0] = org[0] / l;
  org_e[1] = org[1] / l;
  org_e[2] = org[2] / l;
}
 
 double inner(double x[3], double y[3])
{
  return x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
}

void save_file(FILE *f, void *buffer)
{
  printf("%d\n", time);
  /* write PPM file */
  if (f) {
    int i, x, y;
    GLubyte *ptr;
    ptr = buffer;

    fprintf(f,"P6\n");
    fprintf(f,"# ppm-file created\n");
    fprintf(f,"%i %i\n", Width,Height);
    fprintf(f,"255\n");
    for (y=Height-1; y>=0; y--) {
      for (x=0; x<Width; x++) {
        i = (y*Width + x) * 4;
        fputc(ptr[i], f);   /* write red */
        fputc(ptr[i+1], f); /* write green */
        fputc(ptr[i+2], f); /* write blue */
      }
    }
    fclose(f);
  }
  printf("all done\n");
}
