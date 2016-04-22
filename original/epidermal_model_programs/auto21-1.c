/* 
   ATP and B are calculated at each grid point (j, k, l), where x = j * (LX / NX) and so on.
   B.C. are a[NX][k][l]=a[0][k][l], a[j][NY][l]=a[j][0][l], and so on. 
*/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "param.h"
#include "funcs.h"
#include "main.h"

#define EPS 1.0e-7

void error_msg(char ss[]);
void input_uvp(char str[], double u[], double v[], double p[], double a[NX+1][NY+1][NZ+1], 
	       double B[NX+1][NY+1][NZ+1], double w[][N2], int ncell_var);

int main(int argc, char *argv[])
{
	/*
		
	*/
  double u[NN], v[NN], p[NN], fat[NN], Vc[NN], w[NN][N2], a[NX+1][NY+1][NZ+1], B[NX+1][NY+1][NZ+1], u_ave[NN]; //too large stack allocation
  double xx[NN], yy[NN], zz[NN], r[NN];
  double ageb[NN], agek[NN];
  int state[NN], div_times[NN], touch[NN];
  int lj[NN][N2];
  int indx[NN]; // lj[i][0...indx-1]
  int ljd[NN][N2];
  int ljd_indx[NN]; // ljd[i][0...ljd_indx-1]

  /* debug */
  int pair[NN], pair2[NN], pair_indx = 0;
  double L[NN]; 
  int other_cell[NN];
  int tb[NN];
  int ncell_var, num;
  int i, j, k;
  int KEY_INPUT_DATA;
  char celldata_name[100];
  char str[100], str_chem[100];
  FILE *finput, *fdebug;
  void checkstate(int [], int );
  void checkparam();

  int nmx = (int)(COMPRESS_FACTOR * LX/(2.0*R_memb));
  int nmy = (int)(COMPRESS_FACTOR * LY/(2.0*R_memb));

  printf("mkdir %s\n", OUTPUTDIR );
  sprintf(str, "mkdir -p %s", OUTPUTDIR);
  system(str);

  checkparam();

  if (fabs(dx-dh)>EPS || fabs(dy-dh)>EPS || fabs(dz-dh)>EPS) {
    printf("error: dx=%f dy=%f dz=%f.\n", dx, dy, dz);
    printf("grid size must be dh=%f.\n", dh);
    exit(1);
  }
  if (argc <= 1) { 
    printf("input file name required\n");
    exit(1);
  }

  strcpy(celldata_name, argv[1]);

  if((finput = fopen(celldata_name, "r")) == NULL) {
    printf("input file error2\n");
    exit(1);
  }

  if ((num = count_line(finput)) != NN) {
    printf(" error: input data file not consistent.\n");
    printf("input file lines = %d: NN = %d\n", num, NN);
    exit(5);
  }

  fclose(finput);
  if((finput = fopen(celldata_name, "r")) == NULL) {
    printf("input file error2\n");
    exit(1);
  }
  
  initialize_state(state, ageb, agek, fat, Vc);
  /* input vales, initialize gj, return first blank cell's index*/
  initialize_vars(finput, xx, yy, zz, lj, state, r, ageb, agek, div_times, fat, Vc, touch, L, other_cell, tb, &ncell_var, &NDER, &NMEMB);
  printf("variables initialized: NDER=%d NMEMB=%d ncell=%d NN=%d\n", NDER, NMEMB, ncell_var, NN);

  if (nmx != NMX || nmy != NMY) {
    printf("error: number of membrane particles inconsistent with parameter file\n");
    exit(1);
  }

  checkstate(state, ncell_var);
  
  connect_lj(lj, state, ncell_var, xx, yy, zz, r, indx);
  printf("lj connection initialized\n");  

  check_pair(ncell_var, xx, yy, zz, pair, pair2, &pair_indx, L, other_cell);
  
  fclose(finput);
  
  initial_u(u, v, p, a, B, w); // initialization

  KEY_INPUT_DATA = atoi(argv[2]);
  if (KEY_INPUT_DATA != 1 && KEY_INPUT_DATA != 0) {
    printf("error: 2nd argument must be 0 or 1\n");
    exit(3);
  }
  
  if (KEY_INPUT_DATA) {
    printf("chemical data read from recycled_data\n");
    
    sprintf(str_chem, "recycled_data");
    input_uvp(str_chem, u, v, p, a, B, w, ncell_var);
  }

  if (SYSTEM == WHOLE) {
    printf("computing the whole epidermis.\n");
    //    SYSTEM = WHOLE;
  }
  else if (SYSTEM == BASAL) {
    printf("computing only the basal layer and the dermis.\n");
    //    SYSTEM = BASAL;
  }
  else {
    printf("parameter SYSTEM must be 'WHOLE' or 'BASAL'\n");
    exit(1);
  }

  if (KEY_FORCED_SC) 
    printf("forced cornification enabled\n");
  else
    printf("forced cornification disabled\n");

  if (KEY_INPUT_DATA && KEY_FORCED_SC) {
    printf("error: forced cornification must be disabled\n");
    exit(1);
  }
  else if (!KEY_INPUT_DATA && !KEY_FORCED_SC) {
    printf("WARNING: sc formation would take longer without forced cornification.\n");
  }

  if (KEY_DERMAL_CHANGE)
    printf("computing dermal change\n");
  else 
    printf("fixed dermal shape\n");

  printf("divmax=%d, accel_div=%f MALIGNANT=%d\n", 
	 div_max, accel_div, MALIGNANT);
    printf("K_TOTAL=%f, K_DESMOSOME_RATIO=%f\n", K_TOTAL, K_DESMOSOME_RATIO);
  
  evolution(u, v, p, w, a, B, xx, yy, zz, r, ageb, agek, state, div_times, fat, 
	    Vc, touch, lj, indx, &ncell_var, u_ave, pair, pair2,  &pair_indx,  L, other_cell,
            ljd, ljd_indx, tb);  

  printf("finished\n");

  return 0;
}

void checkparam()
{
  double tdiv, tseparate;

  // division period must be greater than the time for the completion of cell division
  tdiv = para_agki_max * (1.0-stoch_div_time_ratio) / (eps_kb * u0);
  tseparate = 2*R_max / eps_L;

  printf("tdiv=%f; tseparate=%f\n", tdiv, tseparate);
  if (tdiv < tseparate) {
    printf("error: division period must be greater than the time for the completion of cell division\n");
    exit(1);
  }
}
void checkstate(int state[], int ncell)
{
  int i;
  printf("ncell=%d\n", ncell);
  for (i=0; i<ncell; i++) 
    if (state[i] == BLANK) {

      printf("input data error\n");
      exit(1);
    }
  if (state[ncell] != BLANK) {
      printf("input data error\n");
      exit(2);
  }
}

void input_uvp(char str[], double u[], double v[], double p[], double a[NX+1][NY+1][NZ+1], 
	       double B[NX+1][NY+1][NZ+1], double w[][N2], int ncell_var)
{
  int i, j, k, n, nnn, dum=-1;
  char ss[256];
  char funcnam[256]="input_uvp";
  FILE *fin_uvp;

  if((fin_uvp = fopen(str, "r")) == NULL) {
    printf("error: cannot open %s\n", str);
    exit(1);
  }
  n = count_line(fin_uvp);
  fclose(fin_uvp);
  fin_uvp=fopen(str, "r");
  if (n != ncell_var + 2*(NX+1)*(NY+1)*(NZ+1) + ncell_var) {
    printf(" error: input data file initial_data is inconsistent\n");
    printf("n=%d ncell_var=%d, sum=%d\n", n, ncell_var, ncell_var + 2*(NX+1)*(NY+1)*(NZ+1) + ncell_var);
    exit(1);
  }
  
  for(i = 0; i < ncell_var; i++) {
    fgets(ss, sizeof(ss), fin_uvp); 
    nnn=sscanf(ss, "%*d %lf %lf %lf", &(u[i]), &(v[i]), &(p[i]));
    if (nnn != 3)
      error_msg(funcnam);
  }

  for(i = 0; i <= NX; i++) 
    for(j = 0; j <= NY; j++) 
      for(k = 0; k <= NZ; k++) {
	fgets(ss, sizeof(ss), fin_uvp); 
	if(sscanf(ss, "%*f %*f %*f %lf", &(a[i][j][k])) != 1) 
	  error_msg("2");
      }

  for(i = 0; i <= NX; i++) 
    for(j = 0; j <= NY; j++) 
      for(k = 0; k <= NZ; k++)  {
	fgets(ss, sizeof(ss), fin_uvp); 
	if(sscanf(ss, "%*f %*f %*f %lf", &(B[i][j][k])) != 1) 
	  error_msg("3");
      }

  for(i = 0; i < ncell_var; i++) {
    fscanf(fin_uvp, "%d", &dum);
    if (dum != i) {
      error_msg("4");
    }
    for (j = 0; j < N2; j++) {
      if(fscanf(fin_uvp, "%lf", &(w[i][j])) != 1) {
	error_msg("5");
      }
    }
  }
  
  fclose(fin_uvp);
}

void error_msg(char ss[])
{

  printf("error in function %s\n", ss);
  exit(1);
}



