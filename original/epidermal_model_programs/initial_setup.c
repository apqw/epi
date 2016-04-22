#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include "param.h"
#include "SFMT.h"
#include "funcs.h"
#include "main.h"

void relaxation(int state[], double rad[], double xx[], double yy[], double zz[], 
		int nfix, int nlast, double Z_LEVEL);

int main(int argc, char *argv[]) {
  int flg, nlast, DER_LAYER;
  int j, k, l, mm, nx, ny, nfix, mm0;
  double diam;
  FILE *fout;
  char str[256];
  double Z_LEVEL;
  double rad[NN], xx[NN], yy[NN], zz[NN], avzz=0.0, avzz_dum;
  int state[NN];
  double pratio = 2*P_MEMB;
  int tb = 0, tb_fix = 0;

  diam = 2*R_memb;

  if (argc != 3) {
    printf("usage: number of FIX, number of dermal layers\n");
    exit(1);
  }
  
  nfix = atoi(argv[1]);
  DER_LAYER = atoi(argv[2]);

  nx = (int)(LX/(R_memb*pratio));
  ny = (int)(LY/(R_memb*pratio));
  printf("nx=%d ny=%d\n", nx, ny);
  
  Z_LEVEL = DER_LAYER * 2.0*R_der;
  
  // membrane
  mm=0;
  for (k=0; k<ny; k++) {
    for (j=0; j<nx; j++) {
      state[mm] = MEMB;
      rad[mm] = R_memb;
      xx[mm] = 0.5*R_memb + j*R_memb*pratio;
      yy[mm] = 0.5*R_memb + k*R_memb*pratio;
      zz[mm] = Z_LEVEL + R_memb;
      mm++;
    }
  }
  printf("membrane cells=%d\n", mm);
  mm0=mm;
  
  // der
  nx = (int)(LX/(2*R_der));
  ny = (int)(LY/(2*R_der));
  for (l=0; l<DER_LAYER; l++) {
    for (k=0; k<ny; k++) {
      for (j=0; j<nx; j++) {
	state[mm] = DER;
	rad[mm] = R_der;
	xx[mm] = R_der + 2*j*R_der;
	yy[mm] = R_der + 2*k*R_der;
	zz[mm] = R_der + 2*l*R_der;
	mm++;
      }
    }
  }
  printf("DER cells=%d\n", mm-mm0);
  
  nlast = mm+nfix;
  
  for (l=mm; l < nlast; l++) {
    state[l] = FIX;
    rad[l] = R_max;
    xx[l] = LX/2.0 + (LX/4.0)*cos(2*M_PI*(l-mm)/nfix);
    yy[l] = LX/2.0 + (LX/4.0)*sin(2*M_PI*(l-mm)/nfix);
    zz[l] = R_max+R_memb+Z_LEVEL+R_memb;
  }

  // relaxation process:
  printf("place stem cells...");
  fflush(stdout);
  relaxation(state, rad, xx, yy, zz, nfix, nlast, Z_LEVEL);
  printf("done.\n");
  
  sprintf(str, "input_%dlayer_nfix%d.dat", DER_LAYER, nfix);
  
  if ((fout=fopen(str, "w")) == NULL) {
    printf("error\n");
    exit(1);
  }
  
  tb_fix = 0;
  for (l=0; l<nlast; l++) {
    if (state[l] == FIX) {
      tb = tb_fix;
      tb_fix++;
    }
    else tb = 0;
    fprintf(fout, "%d %d %f %f %f %f %f %f %f %f %d %f %f %d %f %d %d\n", 
	    l, state[l], rad[l], 0.0, 0.0, 0.0, xx[l], yy[l], zz[l], 0.0, 0, 0.0, 0.0, 0, 0.0, -1, tb);
  }
  for (l=nlast; l<NN; l++) 
    fprintf(fout, "%d %d %f %f %f %f %f %f %f %f %d %f %f %d %f %d %d\n", 
	    l, BLANK, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0, 0.0, -1, 0);
  
  fclose(fout);
}

void relaxation(int state[], double rad[], double xx[], double yy[], double zz[], 
		int nfix, int nlast, double Z_LEVEL)
{

  double dumxx, dumyy, dumzz;
  int j, l;
  double distlj, dum_lj6, LJ1, dum_ljmain, distmin;

  for (j=nlast-nfix; j<nlast; j++) {
    if (state[j] != FIX) {
      printf("error: state must be FIX\n");
      exit(1);
    }

    distmin=10000.0;
    while(1) {
      for (l=0; l<nlast-nfix; l++) {
	if (state[l] != MEMB) continue;
	distlj = dist3(xx[j], yy[j], zz[j], xx[l], yy[l], zz[l]);
	if (distlj < distmin) distmin = distlj;
      }

      if (distmin < 1.1*(rad[j]+rad[l])) break;
      else {
	if (zz[j] < Z_LEVEL-R_max) {
	  printf("zz too low: zz=%f zz-zlevel=%f\n", zz[j], zz[j]-Z_LEVEL);
	exit(1);
	}
	zz[j] -= 0.001;
      }
    }
  }
}

double dist3(double x1, double y1, double z1, double x2, double y2, double z2)
{
  if(PERIODIC == 1) {
    x2 = ret_x(x1, x2);
    y2 = ret_y(y1, y2);
  }
  return norm(x1 - x2, y1 - y2, z1 - z2);
}
double ret_x(double x1, double x2)
/* choose one of LX-periodic values of x2 such that
   |x1-x2| takes minimum
*/
{
  double min, tmp;
  min = fabs(x2 - x1);  
  tmp = x2;
  if(min > fabs(x2 - LX - x1)) {
    min = fabs(x2 - LX - x1);
    tmp = x2 - LX;
  }
  if(min > fabs(x2 + LX - x1)) {
    min = fabs(x2 + LX - x1);
    tmp = x2 + LX;
  }
  return tmp;
}

double ret_y(double x1, double x2)
{
  double min, tmp;
  min = fabs(x2 - x1);
  tmp = x2;
  if(min > fabs(x2 - LY - x1)) {
    min = fabs(x2 - LY - x1);
    tmp = x2 - LY;
  }
  if(min > fabs(x2 + LY - x1)) {
    min = fabs(x2 + LY - x1);
    tmp = x2 + LY;
  }
  return tmp;
}
double norm(double x, double y, double z)
{
  return sqrt(x * x + y * y + z * z);
}
