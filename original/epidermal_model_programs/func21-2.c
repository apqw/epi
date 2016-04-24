#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "param.h"
#include "funcs.h"

#define SQR_RAD 1.4142135623731 // sqrt of 2

void initial_u(double u[], double v[], double p[], double a[NX+1][NY+1][NZ+1], 
	       double B[NX+1][NY+1][NZ+1], double w[NN][N2])

{
  int i, j, k;

  for(i = 0; i < NN; i++) {
    u[i] = u0;
    v[i] = v0;
    p[i] = p0;
    for(j = 0; j < N2; j++) {
      w[i][j] = w0;
    }
  }
  for(i = 0; i <= NX; i++) {
    for(j = 0; j <= NY; j++) {
      for(k = 0; k <= NZ; k++) {
	a[i][j][k] = a0;
	B[i][j][k] = B0;
      }
    }
  }  
}

void save_start(FILE **fp, char *file_name)
{
  if((*fp = fopen(file_name, "w")) == NULL) {
    printf("output file error1\n");
    exit(1);
  }
  fclose(*fp);
  if((*fp = fopen(file_name, "a")) == NULL) {
    printf("output file error2\n");
    exit(1);
  }
}
void save_start_a(FILE **fp, char *file_name)
{
  if((*fp = fopen(file_name, "a")) == NULL) {
    printf("output file error2\n");
    exit(1);
  }
}

void save_end(FILE **fp)
{
  fclose(*fp);
}

/*debug*/
/*細胞分裂中の細胞がいるかを調べる*/
void check_pair(int ncell, double xx[], double yy[], double zz[], 
		int pair[], int pair2[], int *pair_indx, double L[], int other_cell[])
{
  int i, j, k, l;
  *pair_indx = 0;
  for (j = 0; j < NN; j++) {
    pair[j] =  -1;
    pair2[j] =  -1;
    //    L[j] = 0.0;
  }
  for (j = NMEMB+NDER; j < ncell; j++) {
    if(other_cell[j] != -1){
      l = other_cell[j];
      // pairへの格納において重複をさけるためj < l のみ
      if(j < l){
	//        L[j] = dist3(xx[j], yy[j], zz[j], xx[l], yy[l], zz[l]);
	//        L[l] = L[j];
        pair[*pair_indx] = j; 
        pair2[*pair_indx] = l; 
        (*pair_indx) ++ ;
      }
    }
  }
}

void connect_lj(int lj[][N2], int state[], int ncell, double xx[], double yy[], double zz[], 
		double r[], int indx[])
{
  int i, j, k, l, m, o;
  double distlj, eps=1.0e-5;
  int anx, any, anz;
  int aix, aiy, aiz;
  int p, ii, jj, kk, jl, jr, jb, ju;
  static int area[ANX][ANY][ANZ][N3];   
  static int aindx[ANX][ANY][ANZ];
  static int indx_bak[NN];
  //  static int aindx_der[ANX][ANY][ANZ];
  static int flg=0;
  FILE *fdebug;


  /* if (!flg) { */
  /*   check_area_der(state, xx, yy, zz, area, aindx, aindx_der, ncell); */
  /*   flg = 1; */
  /* } */
  
  check_area(state, xx, yy, zz, area, aindx, ncell);
  
  /* connect MEMB only once : no area information needed */
  //周囲の4つと
  for(j=0; j < NMEMB; j++){
    if(flg) continue; 
    indx[j] = 0;

    jj = j % NMX;
    kk = j / NMX;
    if (jj == 0) jl = j + NMX-1;
    else jl = j-1;
    
    if (jj == NMX-1) jr = j - (NMX-1);
    else jr = j+1;
    
    if (kk == 0) jb = j + NMX*NMY - NMX;
    else jb = j-NMX;
    
    if (kk == NMY-1) ju = j - (NMX*NMY - NMX);
    else ju = j + NMX; 

    lj[j][indx[j]] = jl; indx[j]++;
    lj[j][indx[j]] = jr; indx[j]++;
    lj[j][indx[j]] = jb; indx[j]++;
    lj[j][indx[j]] = ju; indx[j]++;

    indx_bak[j] = indx[j];
  }
  
  for(i=0; i < ncell; i++){ // interaction btw MEMB and stems should be counted
    if(state[i] == DISA) continue;

    indx[i] = (state[i] == MEMB ? indx_bak[i] : 0);

    anx = (int) (xx[i]/AREA_GRID);
    any = (int) (yy[i]/AREA_GRID);
    anz = (int) (zz[i]/AREA_GRID);
    
    if (anx >= ANX || any >= ANY || anz >= ANZ || anx < 0 || any < 0 || anz < 0) {
      printf("error in function connect_lj\n");
      printf("i=%d anx=%d any=%d anz=%d xx=%f yy=%f zz=%f\n", i, anx, any, anz, xx[i], yy[i], zz[i]);
      exit(1);
    }
    
    /* if (state[i] == FIX || state[i] == MUSUME) ii = 3; */
    /* else ii = 2; */
    ii = 2;
    for(j = anx - ii; j <= anx + ii; j++) {
      aix = (j + ANX) % ANX;
      for(k = any - ii; k <= any + ii; k++) {
	aiy = (k + ANY) % ANY;
	for(l = anz - ii; l <= anz + ii; l++) {
	  aiz = (l + ANZ) % ANZ;
	  for(m = 0; m < aindx[aix][aiy][aiz]; m++) {
	    o = area[aix][aiy][aiz][m];
	    if((i != o) && state[o] != DISA) {
	      if (state[i] == MEMB && state[o] == MEMB) { // already treated above
		continue;
	      }
	      distlj = dist3(xx[i], yy[i], zz[i], xx[o], yy[o], zz[o]);

	      if(distlj <= LJ_THRESH*(r[i] + r[o])) {
		lj[i][indx[i]] = o;
		indx[i] += 1;
		if (indx[i]>=N2) {
		  printf("error: i=%d state[i]=%d indx[i]=%d: indx must be < %d\n", i, state[i], indx[i], N2);
		  if ((fdebug=fopen("debug_info", "w")) == NULL) {
		    printf("error: cannot open debug_info\n");
		    exit(1);
		  }
		  for (jj=0; jj<indx[i]; jj++)
		    fprintf(fdebug, "%d %d %d\n", jj, lj[i][jj], state[lj[i][jj]]);
		  fclose(fdebug);
		  exit(1);
		}
	      }
	    }
	  }
	}
      }
    }
  }
}

void input_conc(double u[], double v[], double p[], double w[][N2], double a[NX+1][NY+1][NZ+1], 
		double B[NX+1][NY+1][NZ+1], double ini_u[], double ini_v[], double ini_p[], 
		double ini_w[][N2], double ini_a[NX+1][NY+1][NZ+1], double ini_B[NX+1][NY+1][NZ+1])
{
  int i, j, k;
  for(i = 0; i < NN; i++) {
    u[i] = ini_u[i];
    v[i] = ini_v[i];
    p[i] = ini_p[i];
    for(j = 0; j < N2; j++) {
      w[i][j] = ini_w[i][j];
    }
  }
  
  for(i = 0; i <= NX; i++) {
    for(j = 0; j <= NY; j++) {
      for(k = 0; k <= NZ; k++) {
	a[i][j][k] = ini_a[i][j][k];
	B[i][j][k] = ini_B[i][j][k];
      }
    }
  }
}

void check_zzmax(double zz[], int state[], int ncell, double *zzmax)
{
  int j;
  *zzmax = 0.0;
  for (j=NMEMB+NDER;j<ncell; j++) {
    if (state[j] == DEAD || state[j] == ALIVE || state[j] == MUSUME || state[j] == FIX) {
      if (*zzmax < zz[j]) *zzmax = zz[j];
    }
  }
}

/* 基底層の配置と結合を入力 */
/* 細胞の配置と結合を入力 */
void initialize_vars(FILE *fi, double xx[], double yy[], double zz[], 
		     int lj[NN][N2], int state[], double r[], double ageb[], double agek[], int div_times[], 
		     double fat[], double Vc[], int touch[], double L[],
		     int other_cell[], int tb[], int *ncell_var, int *nder, int *nmemb)
{
  int i, j, k;
  int n = 0;
  int n2 = N2;
  int nn = NN;
  int flgder=0;
  int flgair=0;
  
  *nder = 0;
  *nmemb = 0;
  for(i = 0; i < nn; i++) {
    fscanf(fi, "%*d %d %lf %lf %lf %*f %lf %lf %lf %*f %d %lf %lf %d %lf %d %d", 
	   &state[i], &r[i], &ageb[i], &agek[i], &xx[i], &yy[i], &zz[i], 
	   &div_times[i], &fat[i], &Vc[i], &touch[i], &L[i], &other_cell[i], &tb[i]);
    if (SYSTEM == BASAL && (state[i] == ALIVE || state[i] == DEAD || state[i] == AIR)) {
      printf(" input date must not contain ALIVE or DEAD in case of BASAL\n");
      exit(1);
    }
    if (state[i] == DER) {
      if (r[i] != R_der) {
	printf("radii of DER not consistent with param.h\n");
	exit(1);
      }
      (*nder)++;
    }
    else if (state[i] == MEMB) {
      if (r[i] != R_memb) {
	printf("radii of DER not consistent with param.h\n");
	exit(1);
      }
      (*nmemb)++;
    }
    if(n == 0 && state[i] == BLANK) {
      n = i;
    }
  }
  
  for(j = 0 ; j < NN ; j++) {
    for(k = 0 ; k < N2  ; k++)
      lj[j][k] = -1;
  }
  *ncell_var = n;
}

void reopen_file(FILE **f, char *str)
{
  fclose(*f);
  if((*f = fopen(str, "r")) == NULL) {
    printf("input file error2\n");
    exit(1);
  }
}

int count_line(FILE *f)
{
  int c;
  int cnt=0;
  while((c=getc(f)) != EOF) {
    if (c=='\n') cnt++;
  }
  return cnt;
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
/* ノルムを計算 */
double norm(double x, double y, double z)
{
  return sqrt(x * x + y * y + z * z);
}

/* 距離を計算 (2次元) */
double dist(double x1, double y1, double x2, double y2)
{
  if(PERIODIC == 1) {
    x2 = ret_x(x1, x2);
    y2 = ret_y(y1, y2);
  }
  return norm(x1 - x2, y1 - y2, 0.);
}

/* 距離を計算 (3次元) */
double dist3(double x1, double y1, double z1, double x2, double y2, double z2)
{
  if(PERIODIC == 1) {
    x2 = ret_x(x1, x2);
    y2 = ret_y(y1, y2);
  }
  return norm(x1 - x2, y1 - y2, z1 - z2);
}

void output_data(int i, double xx[], double yy[], double zz[], int lj[NN][N2],int state[], double r[], 
		 double ageb[], double agek[], double u[], double w[NN][N2], double av[], 
		 int div_times[], double L[], double fat[], double Vc[], int touch[], int other_cell[], int tb[], int indx[], int ncell, double a[NX+1][NY+1][NZ+1], 
		 double B[NX+1][NY+1][NZ+1])
{
  int j, k, l, m;
  int nn = NN;
  int n2 = N2;

  FILE *fo;
  char filename[256];

  sprintf(filename, "%s/%d", OUTPUTDIR, i);

  if((fo = fopen(filename, "w")) == NULL) {
    printf("output file error3\n");
    exit(1);
  }

  for(j = 0; j < nn; j++) {
    fprintf(fo, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf "
	    "%d %lf %lf %d %lf %d %d\n", 
	    j, state[j], r[j], ageb[j], agek[j], u[j], xx[j], yy[j], zz[j], av[j],  
	    div_times[j], fat[j], Vc[j], touch[j], L[j], other_cell[j], tb[j]);  
  }

  fclose(fo);
}

void initialize_state(int *state, double *ageb, double *agek, double *fat, double *Vc)
{
  int i;
  for(i = 0; i < NN; i++) {
    state[i] = BLANK;
    ageb[i] = 0.;
    agek[i] = 0.;
    fat[i] = 0.0;
    Vc[i] = 0.0;
  }
}

/*areaにどの細胞が乗っているかを記録　   */
/*格子状の点にN3個以上細胞が乗っているときはエラーを返す   */

void check_area(int *state, double *xx, double *yy, double *zz,
		int area[ANX][ANY][ANZ][N3], int aindx[ANX][ANY][ANZ], int ncell)
{
  int i, j, k;
  int aix, aiy, aiz;

  for (i=0; i<ANX; i++)
    for (j=0; j<ANY; j++)
      for (k=0; k<ANZ; k++) aindx[i][j][k] = 0;
  
  for(j = 0; j < ncell; j++) {
    if(state[j] != DISA) {

      aix = (int)(ret_x(LX / 2., xx[j]) / AREA_GRID);
      aiy = (int)(ret_y(LY / 2., yy[j]) / AREA_GRID);
      aiz = (int)((zz[j] < 0. ? 0 : zz[j]) / AREA_GRID);

      if (aix >= ANX || aiy >= ANY || aiz >= ANZ || aix < 0 || aiy < 0 || aiz < 0) {
	printf("error in check_area\n");
	printf("j=%d aix = %d, aiy = %d, aiz = %d, %f %f %f\n",
	       j, aix, aiy, aiz, xx[j], yy[j], zz[j]);
	exit(1);
      }
 
      area[aix][aiy][aiz][aindx[aix][aiy][aiz]] = j;
      aindx[aix][aiy][aiz] += 1;

      if (aindx[aix][aiy][aiz] >= N3) {
	printf("\n");
	printf("j=%d: aix=%d, aiy=%d, aiz=%d, aindx=%d\n",
	       j, aix, aiy, aiz, aindx[aix][aiy][aiz]);
	printf("error in function check_area()\n");
	exit(1);
      }
    }
  }
}

void setup_map(int map[NX+1][NY+1][NZ+1], int map2[NX+1][NY+1][NZ+1], 
	       int state[], int ixx[], int iyy[], int izz[], 
	       double xx[], double yy[], double zz[], double r[], int n)
/*細胞内にある点-> 細胞番号*/
/*細胞外にある点-> -1*/
/* irx = (int) Rad * NX / LENGTH_X = (int) Rad / dx */
{
  int j, k, l, m, ix, iy, iz, ipx, ipy, ipz, imx, imy, imz;
  double mx, my, mz;

  for(j = 0; j <= NX; j++) 
    for(k = 0; k <= NY; k++) 
      for(l = 0; l <= NZ; l++) {
	map[j][k][l] = -1;
	map2[j][k][l] = 0;
      }

  for(j = 0; j < n; j++) {
    if(state[j] != DISA){
      ix = ixx[j];
      iy = iyy[j];
      iz = izz[j];

      // periodic for x and y
      // neumann for z
      // imx, imy, imz... can be out of bound
      // ipx, ipy, ipz... correction for B.C. applied
      for (k = -irx; k <= irx; k++) {
	imx = ix + k;
	mx = imx * dx;
	if (imx < 0) ipx = imx + NX;
	else if (imx >= NX) ipx = imx - NX;
	else ipx = imx;
	for (l = -iry; l <= iry; l++) {
	  imy = iy + l;
	  my = imy * dy;
	  if (imy < 0) ipy = imy + NY;
	  else if (imy >= NY) ipy = imy - NY;
	  else ipy = imy;
	  for (m = -irz; m <= irz; m++) {
	    ipz = imz = iz + m;
	    mz = imz * dz;
	    if ((imz > 0 && imz < NZ) && is_within_cell(mx, my, mz, xx[j], yy[j], zz[j], r[j])) 
	      map[ipx][ipy][ipz] = j;
	  }
	}
      }
      if (state[j] == FIX || state[j] == MUSUME || state[j] == ALIVE) {
	for (k = -2*irx; k <= 2*irx; k++) {
	  imx = ix + k;
	  mx = imx * dx;
	  if (imx < 0) ipx = imx + NX;
	  else if (imx >= NX) ipx = imx - NX;
	  else ipx = imx;
	  for (l = -2*iry; l <= 2*iry; l++) {
	    imy = iy + l;
	    my = imy * dy;
	    if (imy < 0) ipy = imy + NY;
	    else if (imy >= NY) ipy = imy - NY;
	    else ipy = imy;
	    for (m = -2*irz; m <= 2*irz; m++) {
	      ipz = imz = iz + m;
	      mz = imz * dz;
	      if ((imz > 0 && imz < NZ) && 
		is_within_cell(mx, my, mz, xx[j], yy[j], zz[j], FAC_MAP*r[j])) 
		map2[ipx][ipy][ipz] = 1;
	    }
	  }
	}
      }
    }
    
  }

  /* B.C. */
  for (l=0; l<NZ; l++) {
    for (j=0; j<NX; j++) {
      map[j][NY][l] = map[j][0][l];
      map2[j][NY][l] = map2[j][0][l];
    }
    for (k=0; k<=NY; k++) {
      map[NX][k][l] = map[0][k][l];
      map2[NX][k][l] = map2[0][k][l];
    }
  }
}

void output_data_on_map(FILE **fp2, FILE **fp3, int map[NX+1][NY+1][NZ+1], 
			int state[], double u[], double v[], double p[], 
			double agek[], double fat[], 
			double a[NX+1][NY+1][NZ+1],double B[NX+1][NY+1][NZ+1], 
			int kitei_c, int i)
{
  int j, k, m, l;
  double U[NX+1][NY+1][NZ+1], V[NX+1][NY+1][NZ+1], 
    P[NX+1][NY+1][NZ+1],FAT[NX+1][NY+1][NZ+1], AGE[NX+1][NY+1][NZ+1];

  for(j = NX / 2 - irx; j <= NX / 2 + irx; j++) {
    for(k = 0; k <= NY; k++) {
      for(m = 0; m <= NZ; m++) {
	if((l = map[j][k][m]) > -1) {
	  if(state[l] == ALIVE /*|| state[l] == FIX*/) {
	    U[j][k][m] = u[l];
	    V[j][k][m] = v[l];
	    P[j][k][m] = p[l];
	    AGE[j][k][m] = agek[l];
	    FAT[j][k][m] = fat[l];
	  } else if(state[l] == DEAD) {
	    U[j][k][m] = -1.;
	    AGE[j][k][m] = agek[l];
	    FAT[j][k][m] = fat[l];
	  } else if(state[l] == FIX || state[l] == MUSUME) {
	    U[j][k][m] = -2.;
	    AGE[j][k][m] = 0.0;//agek[l];
	    FAT[j][k][m] = 0.0;
	  }
	  else {
	    U[j][k][m] = 0.;
	    AGE[j][k][m] = 0.;
	    FAT[j][k][m] = 0.;
	  }
	} else {
	  U[j][k][m] = 0.;
	  V[j][k][m] = 0.;
	  P[j][k][m] = 0.;
	  AGE[j][k][m] = 0.;
	  FAT[j][k][m] = 0.;
	}
      }
    }
  }
  
  printf("a %.6e  B %.6e  p %.6e  u %.6e  v %.6e DEAD\n", 
	 a[NX / 2][NY / 2][0], B[NX / 2][NY / 2][0], p[kitei_c], u[kitei_c], v[kitei_c]);
  fprintf(*fp2, "%d %.10e %.10e %.10e %.10e %.10e\n", 
	  i, a[NX / 2][NY / 2][0], B[NX / 2][NY / 2][0], p[kitei_c], u[kitei_c], v[kitei_c]);
  fflush(*fp2);
  
  fprintf(*fp3, "\n");
  j = NX / 2;
  k = NY / 2;
  for(l = 0; l <= NZ; l++) 
    fprintf(*fp3, "%d %f %f %f %f %f %f\n", i/CUT, j*dx, 
	    U[j][k][l], V[j][k][l], P[j][k][l], a[j][k][l], B[j][k][l]);
  fflush(*fp3);

}

void convert_coordinate_to_lattice(double xx[], double yy[], double zz[], 
				   int ixx[], int iyy[], int izz[], int ncell)
{
  int j;
  for(j = 0; j < ncell; j++) {
    ixx[j] = ((int)(xx[j] / dx) + NX) % NX;
    iyy[j] = ((int)(yy[j] / dy) + NY) % NY;
    izz[j] = ((int)(zz[j] / dz) + NZ) % NZ;

    if (ixx[j] < 0) { printf("error: invalid index\n"); exit(1); }
    if (iyy[j] < 0) { printf("error: invalid index\n"); exit(1); }
    if (izz[j] < 0) { printf("error: invalid index\n"); exit(1); }
  }
}

void errorcheck(int ncell)
{
  if (ncell >= NN) {
    printf("error: number of cells exceeded the upper bound.\n");
    exit(2);
  }
}

/*
	ncell:セル数？
	nvac
	
	*/
void check_vacant_seats(int ncell, int *nvac, int vanished_cell[], int state[])
{
  int j;
  *nvac = 0; 
  for (j=NDER+NMEMB; j < ncell; j++) {
    if (state[j] == DISA) {
      (*nvac)++;
      vanished_cell[*nvac] = j; // vanished_cell[0] not used
    }
  }
}

#define MAX_COST_FUNC 1.0e2
void calc_energy(int i, double xx[], double yy[], double zz[], int state[], int div_times[],
		 double u[], int touch[], int tb[], int ncell, int nborn, int ndisap, int NFIX)
{
  FILE *fenergy;
  int j, k, l, nlocal=0, nfix=0, nmusume=0, nalive=0, ndead=0, nnodiv=0, nder = 0;
  double zz_av=0.0, u_av=0.0, zz_av2 = 0.0, en_zz=0.0, en_u=0.0, en_zz2 = 0.0;
  double dumzz[ncell], dumu[ncell], dumzz2[ncell];
  //  int nd = 0, ndt = 0;
  int nsub[NFIX];
  
  if ((fenergy=fopen(NAME_ENERGY_OUTPUT, "a")) == NULL) {
    printf("error in calc_energy\n");
    exit(1);
  }
  
  for (j=0; j<NFIX; j++) nsub[j] = 0;

  for (j=NDER+NMEMB; j<ncell; j++) {
    switch(state[j]) {
    case FIX:
      nfix++;
      break;
    case MUSUME:
      nmusume++;
      if (div_times[j] == 0) nnodiv++;
      for (k=0; k<NFIX; k++)  // now the id of FIX is from 0 to nfix-1
	if (tb[j] == k) nsub[k]++;
      break;
    case ALIVE:
      nalive++;
      if (touch[j] == 1) {
	dumzz[nlocal] = zz[j];
	dumu[nlocal] = u[j];
	nlocal++;
      }
      break;
    case DEAD:
      ndead++;
      break;
    default:
      break;
    }
  }
  
  if (nlocal==0) {
    u_av = zz_av = 0.0;
    en_zz = en_u = MAX_COST_FUNC;
  }
  else {
    for(j=0; j < nlocal; j++) {
      zz_av += dumzz[j];
      u_av += dumu[j];
    }
    zz_av /= nlocal;
    u_av /= nlocal;
    
    for(j=0; j < nlocal; j++) {
      en_zz += (dumzz[j]-zz_av)*(dumzz[j]-zz_av);
      en_u += (dumu[j]-u_av)*(dumu[j]-u_av);
    }
    
    en_zz /= (4*Rad*Rad*nlocal);
    en_u /= (u_av*u_av*nlocal);
  }

  if (nfix != NFIX) { 
    printf("error in function calc_energy: NFIX inconsistent\n");
    exit(1);
  } 
    
  fprintf(fenergy, "%f %e %e %e %e %d %d %d %d %d %d %d %d", 
	  i*DT, en_zz, zz_av, en_u, u_av,  
	  nfix, nmusume, nalive, ndead, nlocal, nborn, ndisap, nnodiv);
  for (k=0; k<NFIX; k++) 
    fprintf(fenergy, " %d", nsub[k]);
  fprintf(fenergy, "\n");
  
  
  if (i % (10*CUT) == 0) {
    fprintf(stdout, "time, nfix, nmusume, nalive, nlocal, nborn, ndisap, nnodiv\n");
    fprintf(stdout, "t=%.1f en_zz=%.3e zz_av=%.3e, en_u=%.3e, u_av=%.3e\n"
	    "fix=%d musume=%d alive=%d dead=%d\n"
	    "local=%d born=%d disap=%d nodiv=%d\n"
	    "dmax=%f dmin=%f\n", 
	    i*DT, en_zz, zz_av, en_u, u_av, 
	    nfix, nmusume, nalive, ndead, nlocal, nborn, ndisap, nnodiv, 
	    distmax, distmin);
    fflush(stdout);
  }  
  
  fclose(fenergy);
}

void check_localization(int lj[NN][N2], int indx[], int state[], int touch[], int ncell)
{
  int j, m, l;
  for (j = 0; j<ncell; j++) {
    touch[j] = 0;
    if (state[j] == ALIVE) {
      for(m = 0; m < indx[j]; m++) {
	if ((l=lj[j][m]) > -1) {
          if(state[l]==DEAD) {
            touch[j]=1;
            break;
          }
	}
      }
    }
  }
}

void print_data(int time, int param_div_max, int cut)
{
  printf("total calculation time=%d, cut=%d\n", time, cut);
}

void input_others(double u_ave[], double u[], double v[], double p[], double r[], double xx[], 
		  double yy[], double zz[], double ageb[], double agek[], 
		  int state[], int div_times[], int ncell)
{
  int j;
  for(j = 0; j < NN; j++){
    u_ave[j] = u[j];
  }
  
  for(j = NMEMB + NDER ; j < ncell ; j++){
    if(state[j] == DEAD){
      u[j] = 0.;
      u_ave[j] = 0.;
    }
  }
  
  /*************************************/
  for(j = 0  ; j < NN ; j++){
    if(state[j] == BLANK) {
      u[j] = 0.0;
      v[j] = 0.0;
      p[j] = 0.0;
      u_ave[j] = 0.0;
      r[j] = 0.0;
      xx[j] = 0.0;
      yy[j] = 0.0;
      zz[j] = 0.0;
      
      ageb[j] = 0.0;
      agek[j] = 0.0;

      div_times[j] = 0; 
    }
  }
}

void recycle(int n,double xx[], double yy[], double zz[], int state[], double r[], 
		 double ageb[], double agek[], double av[], 
	     int div_times[], double L[], double fat[],double Vc[],int touch[], 
	     int other_cell[],int tb[],  double u[],double v[],double p[],
	     double w[NN][N2],double a[NX+1][NY+1][NZ+1],double B[NX+1][NY+1][NZ+1])
{
  int j, k, l;
  int nn = NN;
  int n2 = N2;
  int n3 = N3;
  int new_indx[NN];
  int dum_other_cell;

  FILE *fo,*fp1;
  char filename[256];
  sprintf(filename, "%s/%s", OUTPUTDIR, "recycle_cell");
  if((fo = fopen(filename, "w")) == NULL) {
    printf("output file error4\n");
    exit(1);
  }

  k=0;
  for (j = 0; j < nn; j++) {
    if(state[j]!=DISA && state[j]!=BLANK){
      new_indx[j] = k;
      k ++ ;
    }
  }
  k = 0;
  for(j = 0; j < nn; j++) {
	  if(state[j]!=DISA && state[j]!=BLANK){
            if(other_cell[j] !=  -1){
              dum_other_cell = new_indx[other_cell[j]];
            }else{
              dum_other_cell = other_cell[j];
            }
            fprintf(fo, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %d %lf %lf %d %lf %d %d\n", 
		    k, state[j], r[j], ageb[j], agek[j], u[j], xx[j], yy[j], zz[j], av[j], 
		    div_times[j], fat[j], Vc[j], touch[j], L[j], dum_other_cell, tb[j]);
            k++;  
     }
   }
  for(; k < nn; k++) {
    fprintf(fo, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %d %lf %lf %d %lf %d %d\n", 
	    k, BLANK, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0, 0.0, 0,0);
  }
  /*debug*/
  /*もらったプログラムには書いていなかったので書き足した.*/
  save_end(&fo);

  // printf("save start2 ");
  k=0;
  save_start(&fp1, "recycle_last_data_uvp");
  for(l = 0; l < n; l++) {
    if(state[l]!=DISA&&state[l]!=BLANK){
      fprintf(fp1, "%d %.10e %.10e %.10e\n", k++, u[l], v[l], p[l]);
    }
  }
  save_end(&fp1);
  k=0;
  save_start(&fp1, "recycle_last_data_w");
  for(l = 0; l < n; l++) {
    if(state[l]!=DISA&&state[l]!=BLANK){
      fprintf(fp1, "%d ", k++);
      for(j = 0; j < n2; j++) {
	fprintf(fp1, "%lf ", w[l][j]);
      }
      fprintf(fp1, "\n");
    }
  }
  save_end(&fp1);
  save_start(&fp1, "recycle_last_data_a");
  for(l = 0; l <= NX; l++) {
    for(j = 0; j <= NY; j++) {
      for(k = 0; k <= NZ; k++) {
	fprintf(fp1, "%f %f %f %.10e\n", dx*l, dy*j, dz*k, a[l][j][k]);
      }
      fprintf(fp1, "\n");
    }
    fprintf(fp1, "\n");
  }
  save_end(&fp1);
  save_start(&fp1, "recycle_last_data_B");
  for(l = 0; l <= NX; l++) {
    for(j = 0; j <= NY; j++) {
      for(k = 0; k <= NZ; k++) {
	fprintf(fp1, "%f %f %f %.10e\n", dx*l, dy*j, dz*k, B[l][j][k]);
      }
      fprintf(fp1, "\n");
    }
    fprintf(fp1, "\n");
  }
  save_end(&fp1);
}

void initialize_sc(int state[], double zz[], double zzmax, double agek[], int ncell) {
  
  int j;
  
  for (j=NMEMB + NDER; j<ncell; j++) {
    if (state[j] == ALIVE && zzmax - zz[j] < 8 * R_max) {
      agek[j] = (agek[j] > THRESH_DEAD ? agek[j] : THRESH_DEAD);
      state[j] = AIR;
    }
  }
  
  fprintf(stdout, "done.\n");
  fflush(stdout);
}

double positive(double x)
{
  return x > 0.0 ? x : 0.0;
}

int is_within_cell(double x1, double y1, double z1, double x2, double y2, double z2, double r)
/* ecc: degree of deformation (=0 for spheroid) */
{
  // spheroid
  if (dist3(x1, y1, z1, x2, y2, z2) < r) return 1;
  else return 0;

// ellipsoid
  /* double rx = r/obl; */
  /* double ry = r/obl; */
  /* double rz = obl*obl*r; */

  /* if (ellipse_equation(x1, y1, z1, x2, y2, z2, rx, ry, rz) < 1.0) return 1; */
  /* else return 0; */
}

void calc_radius(double *ra, double *rb, double *dum_dist, double x1, double y1, double z1, double r1, 
		 double x2, double y2, double z2, double r2)
{
  double rx1, ry1, rz1, rx2, ry2, rz2;
  void calculate_radius2(double , double , double , double , double , double , 
			 double , double , double , double , double , double , double *, double *);

  *dum_dist = dist3(x1, y1, z1, x2, y2, z2);	  

  // spheroid
  *ra = r1;
  *rb = r2;

  //ellipsoid
  /* rx1 = ry1 = r1/obl1; */
  /* rz1 = r1*obl1*obl1; */
  /* rx2 = ry2 = r2/obl2; */
  /* rz2= r2*obl2*obl2; */

  /* calculate_radius2(x1, y1, z1, x2, y2, z2, rx1, ry1, rz1, rx2, ry2, rz2, ra, rb); */
}

void calculate_radius2(double x1, double y1, double z1, double x2, double y2, double z2, 
		       double rx1, double ry1, double rz1, double rx2, double ry2, double rz2, double *r1, double *r2)
{
  double sinf = 0.0;
  double sin2f = 0.0;
  double cos2f = 0.0;
  double cost = 0.0;
  double sin2t = 0.0;
  double cos2t = 0.0;
  double a = 0.0;
  double b = 0.0;
  sinf = (y2 - ret_y(y2, y1))/dist3(x1, y1, 0.0, x2, y2, 0.0);
  cost = (z2 - z1)/dist3(x1, y1, z1, x2, y2, z2);
  cos2t =cost*cost;
  sin2f =sinf*sinf;
  cos2f = (1.0 - sin2f);
  sin2t = (1.0 - cos2t);
  a = sqrt(sin2t*cos2f/(rx1*rx1) + sin2t*sin2f/(ry1*ry1) + cos2t/(rz1*rz1));
  b = sqrt(sin2t*cos2f/(rx2*rx2) + sin2t*sin2f/(ry2*ry2) + cos2t/(rz2*rz2));
  *r1 = 1.0/a;
  *r2 = 1.0/b;
}

