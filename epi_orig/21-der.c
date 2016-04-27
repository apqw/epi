/*
  NEWモデル [1貯蔵庫バージョン]
  a : ATP 連続場
  p : IP3, u : Ca2+, 
  v : 不活性化変数 離散場[ランダム]

  数値計算 : Euler法
  可視化 : OpenGL
  wの値も出力＆入力する
  aの初期値はファイルから読まない (0で固定)

  Caの計算 + 細胞分裂，移動
  フラットな真皮上での角層形成シミュレーション
  平均化法を用いた数値計算
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "param.h"
#include "funcs.h"



void output2(int num, int ncell, double xx[], double yy[], double zz[], 
	     int state[], double r[], double ageb[], double agek[], double u_ave[], 
	     int div_times[], double fat[], double Vc[], int touch[], double u[], double v[], 
	     double p[], double w[NN][N2], double a[NX+1][NY+1][NZ+1], double B[NX+1][NY+1][NZ+1]);

void evolution(double u[NN], double v[NN], double p[NN], double w[NN][N2], double a[NX+1][NY+1][NZ+1],
	  double B[NX+1][NY+1][NZ+1], double xx[NN], double yy[NN], double zz[NN], 
	  double r[NN], double ageb[NN], double agek[NN], 
	  int state[NN], int div_times[NN], double fat[NN], double Vc[NN], int touch[NN], 
	       int lj[][N2], int indx[], int *ncell_var, double u_ave[NN], int pair[NN], int pair2[NN], int *pair_indx, double L[NN], int other_cell[NN], int ljd[NN][N2], int ljd_indx[], int tb[])
{
  int i, j, k, l;
  double zzmax, debug_cm[3];
  int ncell = *ncell_var, nvac, nborn=0, ndisap=0, vanished_cell[NN];
  // nvac = #vacant seats (present #DISA)
  // ndisap = cumulative #DISA
  
  char str_energy[256];
  int ixx[NN], iyy[NN], izz[NN];
  void dumpall(double [], double [], double [], double []);  
  double old_B[NX+1][NY+1][NZ+1];
  int flg_forced_sc = KEY_FORCED_SC;
  int num_sc = 0;
  int map[NX+1][NY+1][NZ+1], map2[NX+1][NY+1][NZ+1];  
  int sw = INIT_SW;
  int NFIX=0;
  FILE *fdebug_cm; 
  sprintf(str_energy, "mv %s %s.bak", NAME_ENERGY_OUTPUT,  NAME_ENERGY_OUTPUT);
  system(str_energy);

  check_vacant_seats(ncell, &nvac, vanished_cell, state);
  ndisap = nvac;

  input_others(u_ave, u, v, p, r, xx, yy, zz, ageb, agek, state, div_times, ncell);

  check_zzmax(zz, state, ncell, &zzmax);

  print_data(NUM_ITR*DT, div_max, CUT);

  for (j=0; j<NN; j++) {
    if (state[j] == BLANK) break;
    if (state[j] == FIX) {
      //      tb[j] = j;
      NFIX++;
    }
  }

  if ((fdebug_cm=fopen("debug.dat", "w")) == NULL) {
    printf("error: cannot open debug.dat\n");
    exit(1);
  }
  fclose(fdebug_cm);
struct timeval stop,start;

  /* -------------------------------  メインループ ------------------------------------- */
  for(i = 0; i <= NUM_ITR; i++) {
if(i==1)gettimeofday(&start,NULL);
    //        dumpall(xx, yy, zz, u_ave);
    if(i % CUT == 0){
      
      check_localization(lj, indx, state, touch, ncell);
      output_data(i / CUT, xx, yy, zz, lj, state, r, ageb, agek, u, w, u_ave,
		  div_times, L, fat, Vc, touch, other_cell, tb, indx, ncell, a, B);
      
      fprintf(stdout, "%d ", i / CUT); fflush(stdout);
      if (SYSTEM==WHOLE && i % (CUT*10) == 0) {
	output2(i / CUT, ncell, xx, yy, zz, state, r, ageb, agek, u_ave, 
		div_times, fat, Vc, touch, u, v, p, w, a, B);  
      }
    }

    if (i % OUT_ENERGY == 0) {
      check_localization(lj, indx, state, touch, ncell);
      calc_energy(i, xx, yy, zz, state, div_times, u_ave, touch, tb, ncell, nborn, ndisap, NFIX);
    }

    cell_dynamics(state, ageb, agek, xx, yy, zz, r, 
		  lj, indx, u, v, p, fat, Vc, w, u_ave, div_times, &ncell, &nvac, vanished_cell, i, 
		  &sw, &nborn, &ndisap, touch, pair, pair2, pair_indx,  L, other_cell, ljd, ljd_indx, tb);
    gettimeofday(&stop,NULL);
    printf("cell_dyn_per_sec:%lf\n",(double)i/(stop.tv_sec-start.tv_sec));
    check_zzmax(zz, state, ncell, &zzmax);
    
    if (SYSTEM==WHOLE ) {
      convert_coordinate_to_lattice(xx, yy, zz, ixx, iyy, izz, ncell);
      setup_map(map, map2, state, ixx, iyy, izz, xx, yy, zz, r, ncell);
      
      calc_b(old_B, state, agek, map, map2, zzmax, B);
      
      if(i*DT > T_TURNOVER && flg_forced_sc) {    /*- カルシウムダイナミクス (ca dynamics) */
	flg_forced_sc = 0;
	fprintf(stdout, "t=%f: forced cornification...", i*DT);
	fflush(stdout);
	initialize_sc(state, zz, zzmax, agek, ncell);
	num_sc = NUM_SC_INIT;
      }
      
      if(sw >= SW_THRESH || num_sc > 0) {    /*- カルシウムダイナミクス (ca dynamics) */
	fprintf(stdout, "t=%f: calculating ca_dynamics...", i*DT);
	fflush(stdout);
	ca_dynamics(u_ave, B, u, v, p, a, w, ixx, iyy, izz, lj, indx, state, ageb, agek, map,
		    map2, zzmax, r, xx, yy, zz, ncell);
	fprintf(stdout, "done.\n");
	fflush(stdout);
	if (num_sc > 0) num_sc--;
	sw = 0;
      }
    }
    
  }
  *ncell_var = ncell;
  
}

/*-----------------------*/
  
  
void dumpall(double xx[], double yy[], double zz[], double u[])
{
  int j;
  double dum = 0.0;

  for (j=NDER; j<NN; j++) {
    dum += xx[j]*xx[j] + yy[j]*yy[j] + zz[j]*zz[j] + u[j]*u[j];
  }
  printf("dumpall: %.10e\n", dum);
}

void output2(int num, int ncell, double xx[], double yy[], double zz[], 
	     int state[], double r[], double ageb[], double agek[], double u_ave[], 
	     int div_times[], double fat[], double Vc[], int touch[], double u[], double v[], 
	     double p[], double w[NN][N2], double a[NX+1][NY+1][NZ+1], double B[NX+1][NY+1][NZ+1])
{
  FILE *fp1; 
  int i, j, k;
  char str_uvp[256], str_w[256], str_a[256], str_B[256];
  fprintf(stdout, "ca data output..."); 
  
  system("rm -f last_data*");
  sprintf(str_uvp, "last_data_uvp");
  save_start(&fp1, str_uvp);
  for(i = 0; i < ncell; i++) {
    fprintf(fp1, "%d %.10e %.10e %.10e\n", i, u[i], v[i], p[i]);
  }
  save_end(&fp1);
  // printf("save start2 ");
  sprintf(str_w, "last_data_w");
  save_start(&fp1, str_w);

  for(i = 0; i < ncell; i++) {
    fprintf(fp1, "%d ", i);
    for(j = 0; j < N2; j++) {
      fprintf(fp1, "%lf ", w[i][j]);
    }
    fprintf(fp1, "\n");
  }
  save_end(&fp1);
  
  sprintf(str_a, "last_data_a");
  save_start(&fp1, str_a);

  for(i = 0; i <= NX; i++) {
    for(j = 0; j <= NY; j++) {
      for(k = 0; k <= NZ; k++) {
  	fprintf(fp1, "%f %f %f %.10e\n", dx*i, dy*j, dz*k, a[i][j][k]);
      }
      fprintf(fp1, "\n");
    }
    fprintf(fp1, "\n");
  }
  save_end(&fp1);
  
  sprintf(str_B, "last_data_B");
  save_start(&fp1, str_B);
  for(i = 0; i <= NX; i++) {
    for(j = 0; j <= NY; j++) {
      for(k = 0; k <= NZ; k++) {
  	fprintf(fp1, "%f %f %f %.10e\n", dx*i, dy*j, dz*k, B[i][j][k]);
      }
      fprintf(fp1, "\n");
    }
    fprintf(fp1, "\n");
  }
  save_end(&fp1);

  fprintf(stdout, "done.\n");
}

