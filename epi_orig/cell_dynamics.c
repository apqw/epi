#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "SFMT.h"
#include "param.h"
#include "funcs.h"

#define BETA 0.1
#define STOCHASTIC 1
#define Kspring 25.0 // for FIX
#define Kspring_d 5.0 // for MUSUME
#define Kspring_division 5.0 

#define K_DESMOSOME (K_TOTAL * K_DESMOSOME_RATIO)


#define SQR_RAD 1.4142135623731 // sqrt of 2

FILE *fdebug1, *fdebug2;

void cell_dynamics(int *state, double *ageb, double *agek, 
		   double *xx, double *yy, double *zz, double *r,
		   int lj[NN][N2], int indx[], 
		   double *u, double *v, double *p, double *fat, double *Vc, 
		   double w[NN][N2], double *u_ave, int *div_times, 
		   int *ncell, int *nvac, int vanished_cell[], int i, int *sw, 
		   int *nborn, int *ndisap, int *touch, int *pair, int *pair2, 
		   int *pair_indx, double L[NN], int other_cell[NN],int ljd[NN][N2],int ljd_indx[], int tb[])
{    
  
  double dum_lj3, dum_lj4, dum_lj6, dum_ljmain, dum_norm;
  int old_state[NN];
  double old_xx[NN], old_yy[NN], old_zz[NN], old_r[NN];
  double old_ageb[NN], old_agek[NN], old_fat[NN], old_Vc[NN];
  double x1, y1, z1, x2, y2, z2, x4, y4, z4, r1, r2, r4, ra, rb;
  double x_out, y_out, z_out, xx4, yy4, zz4;
  int j, k, l, o, m, p1, p2, p3, nunpaired, jj, kk;
  double LJ1, LJ2, rsum2, tmp, force, nat_lgth, dist_jl, distlj;
  double dum_z, sum_z = 0.0, zsum;
  double spring_force, dumxx, dumyy, dumzz, rder;
  FILE *fp1;
  double fx1, fy1, fz1, fx2, fy2, fz2, fx3, fy3, fz3;
  int t1, t2, t3;
  double dd_const;
  int seed;
  static int flg=0;
  static int flg_vol = 0;
  double derx[NMEMB/2], dery[NMEMB/2], derz[NMEMB/2];
  static double voln_init;
  double vol, FAC;
  double DD_LJ, dmin;
  static int jr[NMX*NMY], jl[NMX*NMY], jll[NMX*NMY], 
    ju[NMX*NMY], jb[NMX*NMY], jbb[NMX*NMY];
  double dn[NMX*NMY], dm[NMX*NMY], nx[NMX*NMY], mx[NMX*NMY], 
    ny[NMX*NMY], my[NMX*NMY], nz[NMX*NMY], mz[NMX*NMY], ipn[NMX*NMY], ipm[NMX*NMY];
  double cr_dist, lambda_dist;

  if (!flg) {
    //   seed=time(NULL);
    seed = 1;
    init_gen_rand(seed);
    flg = 1;

    setup_memb_indices(jr, jl, jll, ju, jb, jbb);
  }  
  initialize_old_values(old_xx, old_yy, old_zz, old_state, old_ageb, old_agek, 
			old_r, old_fat, old_Vc, xx, yy, zz, state, ageb, agek, r, fat, Vc, *ncell);
  
  // bending 

  distmax=0.0;
  distmin=100.0;

  bend_interac(old_xx, old_yy, old_zz, ipn, ipm, 
	       dn, dm, nx, ny, nz, mx, my, mz, jr, jl, jll, ju, jb, jbb, &distmax, &distmin);

  for (j=0; j<NMEMB; j++) {
    xx[j] += DT * bend_force_sqr(nx, mx, ipn, ipm, dn, dm, j, jr[j], jl[j], jll[j], ju[j], jb[j], jbb[j]);
    yy[j] += DT * bend_force_sqr(ny, my, ipn, ipm, dn, dm, j, jr[j], jl[j], jll[j], ju[j], jb[j], jbb[j]);
    zz[j] += DT * bend_force_sqr(nz, mz, ipn, ipm, dn, dm, j, jr[j], jl[j], jll[j], ju[j], jb[j], jbb[j]);
  }

  /* particle interactions */  

  // memb-memb, memb-wall
  for(j = 0; j < NMEMB; j++){ 
    // memb-wall
    if (old_zz[j] < old_r[j])  //interaction with z=0 plane (equiv. mirror image at (x, y, -z))
      interac_wall(old_zz[j], old_r[j], &(zz[j]));
    
    // memb-memb
    for(k = 0; k < indx[j]; k++) {
      l = lj[j][k];
      if (old_state[l] != MEMB || j>l) continue;
      //if (j > l) continue; // avoid double counting
      
      if (k >= N_NEIGHBOR) {
	printf("error: membrane particles should be labeled k < NN\n");
	exit(1);
      }
      interac_memb_memb(j, l, old_xx, old_yy, old_zz, old_r, old_state, xx, yy, zz);
    }
  }    
  
  //der-der, der-memb, der-wall
  for(j = NMEMB; j < NMEMB+NDER; j++) {
    
    // der-wall
    if (old_zz[j] < old_r[j]) 
      interac_wall(old_zz[j], old_r[j], &(zz[j]));
    
    //der-der, der-others
    for(k = 0; k < indx[j]; k++) {
      l = lj[j][k];

      if (old_state[l] == MEMB || old_state[l] == FIX || old_state[l] == MUSUME 
	  || old_state[l] == ALIVE || old_state[l] == DEAD || old_state[l] == AIR) 
	interac_lj(j, l, old_xx, old_yy, old_zz, old_r, xx, yy, zz); // der-others

      else if (old_state[l] == DER) {  // der-der
	if(l > j) // avoid double counting
	  interac_der_der(j, l, old_xx, old_yy, old_zz, old_r, xx, yy, zz);
      }
      else {
	printf("error: state[%d]=%d state[%d]=%d; illegal interaction\n", j, l, state[j], state[l]);
	exit(1);
      }
    }
  }

  // alive-others, stem-membrane, stem-stem: force to the membrance is calculated as a counter force:
  for(j = NMEMB+NDER; j < *ncell ; j++){ 
    switch(old_state[j]) {
    case DISA:
      break;
      
    case ALIVE: case DEAD: case AIR:
      for (k=0; k<indx[j]; k++) {
	l = lj[j][k];
	if ((old_state[l] == ALIVE || old_state[l] == DEAD || old_state[l] == AIR ) && l > j) continue; // avoid double counting
	if (old_state[l] == DER) continue; // already dealt with
	// basal-supra, supra-supra
	interac_supra_others(j, l, old_xx, old_yy, old_zz, old_r, old_agek, xx, yy, zz);
      }
      break;
      
    case FIX: case MUSUME:
      set_spring_force(j, old_state, other_cell, div_times, &spring_force);
      find_dermis(j, xx, yy, zz, old_state, lj, &p1, indx);

      for(k = 0; k < indx[j]; k++) {
	l = lj[j][k];
	switch(old_state[l]) {
	case ALIVE: case DEAD: case AIR:// already dealt with
	  break;

	case MEMB:
	  if (l == p1) 
	    interac_stem_membrane(j, l, old_xx, old_yy, old_zz, old_r, old_state, spring_force, xx, yy, zz);
	  else 
	    interac_lj(j, l, old_xx, old_yy, old_zz, old_r, xx, yy, zz);
	  break;

	case FIX: case MUSUME:
	  if (j < l) // avoid double counting
	    interac_stem_stem(j, l, old_xx, old_yy, old_zz, old_r, other_cell, xx, yy, zz);
	  break;

	case DER: // already dealt with
	  break;

	default:
	  printf("error: state[%d]=%d state[%d]=%d; illegal interaction\n", j, l, state[j], state[l]);
	  exit(1);
	}
      }
      break;
    }
  }
  
  /* interaction between pairs */  
  for (k = 0; k < *pair_indx; k++) { // no double counting occurs here
    j = pair[k];
    l = pair2[k];
    nat_lgth = L[j];
    dist_jl = dist3(old_xx[j], old_yy[j], old_zz[j], old_xx[l], old_yy[l], old_zz[l]);
    
    force = Kspring_division * (dist_jl - nat_lgth); 
    
    dumxx = DT  * force * (old_xx[j] - ret_x(old_xx[j], old_xx[l]) )/ dist_jl;
    dumyy = DT  * force * (old_yy[j] - ret_y(old_yy[j], old_yy[l]) )/ dist_jl; 
    dumzz = DT  * force * (old_zz[j] - old_zz[l])/ dist_jl;;
    xx[j] -= dumxx;
    xx[l] += dumxx;
    yy[j] -= dumyy;
    yy[l] += dumyy;
    zz[j] -= dumzz;
    zz[l] += dumzz;
  }
  
  // for fixed membrane
  if (KEY_DERMAL_CHANGE==0) {
    for (j = 0; j <NMEMB; j++) {
      xx[j] = old_xx[j];
      yy[j] = old_yy[j];
      zz[j] = old_zz[j];
    }
  }
  
  /*--------------------- State renewal -------------------------*/
  
  /*ageb[j]>=para_agki_maxでFIX,MUSUMEは分裂                      */
  /* contact inhibition not considered */
  /****************************************************************/
  
  for(j = *ncell-1 ; j >= NMEMB+NDER; j--) {  
    switch (old_state[j]) {
    case FIX: 
      find_dermis(j, xx, yy, zz, old_state, lj, &p1, indx);
      if (p1 == -1) {
	printf("error: FIX cell cannot find the dermis\n");
	printf("id=%d xx=%f yy=%f zz=%f\n", j, xx[j], yy[j], zz[j]);
	exit(1);
      }
      if (old_ageb[j] >= para_agki_max_fix*(1.0-stoch_div_time_ratio)) {
	divide_cell(j, xx, yy, zz, r, old_state, state, div_times, u, u_ave, v, p, 
		    old_ageb, old_agek, ageb, agek, old_fat, fat, old_Vc, Vc, ncell, nvac, nborn, vanished_cell, 
		    div_max, FIX, lj, indx, pair, pair2, pair_indx, L, other_cell, tb, p1);
      }
      else {
	ageb[j] = old_ageb[j] + DT * ageb_const(old_state[j], u_ave[j], tb[j]);
      }
      break;
      
    case MUSUME:
      find_dermis(j, xx, yy, zz, old_state, lj, &p1, indx);
      if(p1 == -1 && other_cell[j] == -1) { 
	// differentiate if dermal cell is not found AND unpaired
	if (SYSTEM==WHOLE)
	  state[j] = ALIVE;
	else if (SYSTEM==BASAL) {
	  state[j] = DISA;
	  (*nvac)++;
	  vanished_cell[*nvac] = j;
	  (*ndisap)++;
	}
	break;
      }
      if(div_times[j] > 0 && old_ageb[j] >= para_agki_max*(1.0-stoch_div_time_ratio)) { // div_times inherited
	divide_cell(j, xx, yy, zz, r, old_state, state, div_times, u, u_ave, v, p, 
		    old_ageb, old_agek, ageb, agek, old_fat, fat, old_Vc, Vc, ncell, nvac, nborn, vanished_cell, 
		    div_times[j], MUSUME, lj, indx, pair, pair2, pair_indx, L, other_cell, tb, p1);
      }
      else {
	ageb[j] = old_ageb[j] + DT * ageb_const(old_state[j], u_ave[j], tb[j]);
      }
      // check differentiation
      break;

    case DEAD: case AIR:
      if(old_agek[j] >= ADHE_CONST  && indx[j] <= Nc ) { /* adhe() not depend on fat */
	state[j] = DISA;
	(*nvac)++;
	vanished_cell[*nvac] = j; /* store the label of disappered cells */
	
	(*ndisap)++; // count the number of removed cells
	
      } else {
	agek[j] = old_agek[j] + DT * agek_const(old_state[j], u0, tb[j]);
      }
      break;
      
    case ALIVE: 
      /*debug*/
      if(old_agek[j] >= THRESH_DEAD ) {
	state[j] = DEAD;
	*sw += 1;
      } 
      else {
	agek[j] = old_agek[j] + DT * agek_const(old_state[j], u_ave[j], tb[j]);
	tmp = k_lipid_release(u_ave[j], old_agek[j])*old_Vc[j];
	Vc[j] = old_Vc[j] + DT * (k_lipid(u_ave[j], old_agek[j]) * (1.0 - old_Vc[j]) - tmp);
	fat[j] = old_fat[j] + DT * tmp;
      }
      break;
    }
  }
  
  /*debug*/
  /*細胞分裂時における中心間距離の時間変化と分離*/
  /*※ 相互作用の計算ではL[j]しか使わない*/
  nunpaired=0;
  for (k = 0; k < *pair_indx; k++) {
    j = pair[k];
    l = pair2[k];
    if(L[j] < 2.0*old_r[j]){
      L[j] += DT * eps_L;
      L[l] = L[j];
    } 
    else if ((dist_jl = dist3(xx[j], yy[j], zz[j], xx[l], yy[l], zz[l])) > 0.9 * (old_r[j] + old_r[l])) {
      L[j] = 0;
      L[l] = L[j];
      other_cell[j] = -1;
      other_cell[l] = -1;
      pair[k] = -1;
      pair2[k] = -1;
      nunpaired++;
    }
  }
  if (nunpaired > 0) {
    moveover(pair_indx, pair, pair2, nunpaired);
    *pair_indx -= nunpaired;
  }

  periodic_cond(xx, yy, *ncell);
  
  /*  各細胞のLJをつなぎ直す */
  connect_lj(lj, state, *ncell, xx, yy, zz, r, indx);
  
}

void divide_cell(int j, double xx[], double yy[], double zz[], double r[], 
		 int old_state[], int state[], int div_times[], double u[], double u_ave[], 
		 double v[], double p[], double old_ageb[], double old_agek[], 
		 double ageb[], double agek[], double old_fat[], double fat[], double old_Vc[], double Vc[], 
		 int *ncell, int *nvac, int *nborn,
		 int vanished_cell[], int divtimes_inherited, int i, int lj[][N2], int indx[], int *pair, int *pair2, int *pair_indx, double L[NN], int other_cell[], int tb[], int p1)
{

  int nnew;
  double rh, cos_phase, sin_phase, rand_phase;
  double x1, y1, z1, r1;
  double nx, ny, nz; /* normal vector */
  int dbg;
  static int seed;
  static int flg=0;
  double div_beta, div_gamma, par, rand_phi;
  /*debug*/
  /*分離し終わっていない細胞は分裂させない*/

  /* if(other_cell[j] !=  -1){ */
  /*   printf("%d other_cell=%d, L=%f ageb[j]=%f\n", j, other_cell[j], L[j], ageb[j]); */
  /*   exit(1); */
  /* } */
  if (other_cell[j] != -1) return;

  if (i == MUSUME) par = para_agki_max;
  else par = para_agki_max_fix;
  
  if (tb[j] < MALIGNANT)
    div_beta = accel_div * eps_kb*u0 / (par * stoch_div_time_ratio);
  else 
    div_beta = eps_kb*u0 / (par * stoch_div_time_ratio);

  div_gamma = DT*div_beta;
  
  /* division as a Poisson process */

   if (STOCHASTIC && genrand_res53() > div_gamma) { 
     return; 
   } 

  (*nborn)++;

  if (*nvac == 0) {
    nnew = *ncell;
    (*ncell)++;
    errorcheck(*ncell);

  }
  else {
    nnew = vanished_cell[*nvac];
    (*nvac)--;
  }
  
  other_cell[j] = nnew;
  other_cell[nnew] = j;
  /* debug */
  /* pair[]には数字が小さい方を格納する*/
  /* pair2[]には数字が大きい方を格納する*/
  if(j < other_cell[j]){
    pair[*pair_indx] = j;
    pair2[*pair_indx] = other_cell[j];
  }else{
    pair[*pair_indx] = other_cell[j];
    pair2[*pair_indx] = j;
  }
  L[pair[*pair_indx]] = delta_L;
  L[pair2[*pair_indx]] = delta_L;
  (*pair_indx) ++ ;

  x1=xx[j];
  y1=yy[j];
  z1=zz[j];
  /*debug*/
  /*分裂方式の変更に伴い変更*/
  r1=delta_L;
  //  r1=r[j]; 
  rh = 0.5*r1;
 
  r[nnew] = r[j];

  /* r[j] = rh; */
  div_times[nnew] = divtimes_inherited;
  
  if (i==MUSUME) {
    div_times[j]--;
    div_times[nnew]--;
  }
  
  // direction of cell division 
  div_direction(j, xx, yy, zz, r, old_state, indx, &nx, &ny, &nz, p1);  

  xx[j] = x1 + nx*rh;
  yy[j] = y1 + ny*rh;
  zz[j] = z1 + nz*rh;

  xx[nnew] = x1 - nx*rh;
  yy[nnew] = y1 - ny*rh;
  zz[nnew] = z1 - nz*rh;

  set_newstate(&(old_state[nnew]), &(state[nnew]));

  u[nnew] = u[j];
  u_ave[nnew] = u[j];
  v[nnew] = v[j];
  p[nnew] = p[j];
  
  tb[nnew] = tb[j];
  
  ageb[nnew] = ageb[j] = 0.0;
  old_ageb[nnew] = ageb[nnew] = 0.0;

  old_agek[nnew] = agek[nnew] = 0.0;
  old_fat[nnew] = fat[nnew] = 0.0;
  old_Vc[nnew] = Vc[nnew] = 0.0;
}
               
void periodic_cond(double xx[], double yy[], int ncell)
{
  int j;

  for(j=0;j < ncell;j++){
    if(xx[j] > LX){
      xx[j] = xx[j] - LX;
    }
    else if(xx[j] < 0.0){
      xx[j] = xx[j] + LX;
    }

    if(yy[j] > LY){
      yy[j] = yy[j] - LY;
    }
    else if(yy[j] < 0.0){
      yy[j] = yy[j] + LY;
    }
  }
}

 void initialize_old_values(double *old_xx, double *old_yy, double *old_zz, int *old_state, 
    double *old_ageb, double *old_agek, double *old_r, double *old_fat, double *old_Vc, double *xx, double *yy, double *zz, int *state, 
    double *ageb, double *agek, double *r, double *fat, double *Vc, int ncell)
{
  int j;

  for(j = 0; j < ncell; j++) { // DER must be renewed too
    old_xx[j] = xx[j];
    old_yy[j] = yy[j];
    old_zz[j] = zz[j];
    old_state[j] = state[j];
    old_ageb[j] = ageb[j];
    old_agek[j] = agek[j];
    old_r[j] = r[j];
    old_fat[j] = fat[j];
    old_Vc[j] = Vc[j];
  }

}

void check_lj_connection(double xx[], double yy[], double zz[], int lj[][N2],int indx[], int ncell)
{
  int i, j, k, l, m, o;
  for (j = 0; j < ncell; j++) {
    printf(" %d check start ...\n", j);
    for (k = 0; k < indx[j]; k++) {
      l = lj[j][k];
      for (m = 0; m < indx[l]; m++) {
        o = lj[l][m];
        if(o == j){
          break;
        }
        if(m == indx[l]){
          printf("cannnot equal find %d %d\n", j, l);
          exit(1);
        }
      }
    }
  }
}

void swap_int(int *x, int *y)
{
  int tmp;
  tmp = *x;
  *x = *y;
  *y = tmp;
}

 void moveover(int *pair_indx, int pair[], int pair2[], int nunpaired)
 {
  int j, jj, top_indx=0, num = nunpaired;
  void swap_int(int *, int *);
  int tmp[*pair_indx], tmp2[*pair_indx];

  for (j=0, jj=0; j<*pair_indx; j++) {
    if (pair[j] != -1) {
      tmp[jj] = pair[j];
      tmp2[jj] = pair2[j];
      jj++;
    }
  }
  for (j=0; j<*pair_indx-nunpaired; j++) {
    pair[j] = tmp[j];
    pair2[j] = tmp2[j];
  }
  for (j=*pair_indx-nunpaired; j<*pair_indx; j++) {
    pair[j] = -1;
    pair2[j] = -1;
  }

}

double bend_force(double n[], double m[], double ipn[], double ipm[], double dn[], double dm[], 
	   int j, int jr, int jl, int jll, int ju, int jb, int jbb)
{
  return KBEND * (
		  -(n[jr] - ipn[j]*n[j])/dn[j] - (n[jl] - ipn[jl]*n[j])/dn[j]
		  + (n[j] - ipn[jl]*n[jl])/dn[jl] + (n[jll] - ipn[jll]*n[jl])/dn[jl]
		  - (m[ju] - ipm[j]*m[j])/dm[j] - (m[jb] - ipm[jb]*m[j])/dm[j]
		  + (m[j] - ipm[jb]*m[jb])/dm[jb]  + (m[jbb] - ipm[jbb]*m[jb])/dm[jb]);
}

double bend_force_sqr(double n[], double m[], double ipn[], double ipm[], double dn[], double dm[], 
	   int j, int jr, int jl, int jll, int ju, int jb, int jbb)
{
  return KBEND * (
		  -(1.0-ipn[j])*(n[jr] - ipn[j]*n[j])/dn[j] 
		  +(1.0-ipn[jl])*(-(n[jl] - ipn[jl]*n[j])/dn[j]
				  +(n[j] - ipn[jl]*n[jl])/dn[jl])
		  + (1.0-ipn[jll])*(n[jll] - ipn[jll]*n[jl])/dn[jl]
		  
		  -(1.0-ipm[j])*(m[ju] - ipm[j]*m[j])/dm[j] 
		  +(1.0-ipm[jb])*(-(m[jb] - ipm[jb]*m[j])/dm[j]
				  +(m[j] - ipm[jb]*m[jb])/dm[jb])
		  +(1.0-ipm[jbb])*(m[jbb] - ipm[jbb]*m[jb])/dm[jb]);
}

/* return x1-x2 taking into account periodic bc */
double perio_diff(double x1, double x2, double lb) 
{
double diff=x1-x2;
  if (diff*diff < 0.25*lb*lb) return diff;
  else {
    if (x1 > x2) {
      return diff-lb;
    }
    else {
      return diff+lb;
    }
  }
}

void calc_b(double old_B[NX+1][NY+1][NZ+1], int *state, double *agek,     
    int map[NX+1][NY+1][NZ+1], int map2[NX+1][NY+1][NZ+1], 
    double zzmax, double B[NX+1][NY+1][NZ+1])
{
  int j, k, m, l;
  int prev_x, prev_y, prev_z, next_x, next_y, next_z;
  int iz_bound = (int)((zzmax + FAC_MAP*R_max)/dz); 
  int flg_cornified;
  double dum_age;

  for(j = 0 ; j <= NX ; j++){
    for(k = 0 ; k <= NY ; k++){
      for(l = 0 ; l <= iz_bound; l++){
        old_B[j][k][l] = B[j][k][l];
      }
    }
  }

  for(j = 0; j < NX; j++) {
    prev_x = (j - 1 + NX) % NX;
    next_x = (j + 1 + NX) % NX;
    for(k = 0; k < NY; k++) {
      prev_y = (k - 1 + NY) % NY;
      next_y = (k + 1 + NY) % NY;
      for(l = 0; l < iz_bound; l++) {
        if(l == 0) {
          prev_z = 1;
        } else {
          prev_z = l - 1;
        }
        if(l == NZ) {
          next_z = NZ - 1;
        } else {
          next_z = l + 1;
        }	  
        if ((m = map[j][k][l]) > -1) {
          if (state[m] == ALIVE || state[m] == DEAD) { 
            dum_age = agek[m];
            flg_cornified = 1;
          }
        }
        else {
          dum_age = 0.0;
          flg_cornified = 0;
        }

        B[j][k][l] = old_B[j][k][l] 
          + DT * (DB * (map2[prev_x][k][l] * (old_B[prev_x][k][l] - old_B[j][k][l])
                + map2[j][prev_y][l] * (old_B[j][prev_y][l] - old_B[j][k][l])
                + map2[j][k][prev_z] * (old_B[j][k][prev_z] - old_B[j][k][l])
                + map2[next_x][k][l] * (- old_B[j][k][l] + old_B[next_x][k][l])
                + map2[j][next_y][l] * (- old_B[j][k][l] + old_B[j][next_y][l])
                + map2[j][k][next_z] * (- old_B[j][k][l] + old_B[j][k][next_z])) / (dz * dz)
              + fB(dum_age, old_B[j][k][l], flg_cornified));
      }
    }
  }
  /* B.C. */
  for (l=0; l<=iz_bound; l++) {
    for (j=0; j < NX; j++) B[j][NY][l] = B[j][0][l];
    for (k=0; k <= NY; k++) B[NX][k][l] = B[0][k][l];
  }

}

/* reaction terms for B */
double fB(double age, double B, int flg_cornified)
{
  return IR(age, flg_cornified) - para_kb*B;
}

double IR(double age, int is_cornified) //角化前後で刺激物質を出す
{
  if(is_cornified) {
    if(age > THRESH_DEAD - DUR_ALIVE  && age <= THRESH_DEAD + DUR_DEAD )
      return 1.0;
    else return 0.0;
  }
  else {
    return 0.0;
  }
}

void set_spring_force(int j, int old_state[], int other_cell[], int div_times[], double *spring_force)
{
  if (old_state[j] == FIX || (old_state[j] == MUSUME && other_cell[j] != -1 
			      && old_state[other_cell[j]] == FIX))
    *spring_force = Kspring;
    else if (div_times[j] > 0)
      *spring_force = Kspring_d;
    else *spring_force = 0.0;
}

void interac_der_der(int j, int l, double old_xx[], double old_yy[], double old_zz[], double old_r[], 
		     double xx[], double yy[], double zz[])
{
  double distlj, LJ1, dum_lj6, dum_ljmain, dumxx, dumyy, dumzz;
  distlj = dist3(old_xx[j], old_yy[j], old_zz[j], old_xx[l], old_yy[l], old_zz[l]);
  
  if(distlj < old_r[j] + old_r[l] - delta_R){
    LJ1 = (old_r[j]+old_r[l])/(distlj + delta_R);
    dum_lj6 = pow(LJ1, 6);
    dum_ljmain = 4.0*eps_m * dum_lj6 * (dum_lj6-1.0) /((distlj + delta_R)*distlj) + para_ljp2;
  }
  else if(distlj < old_r[j] + old_r[l]){
    dum_ljmain = para_ljp2;
  }
  else{
    LJ1 = (old_r[j]+old_r[l])/(distlj + delta_R);
    dum_lj6 = pow(LJ1, 6);
    dum_ljmain = 4.0*eps_m * dum_lj6 * (dum_lj6-1.0) / (distlj*distlj) + para_ljp2;
  }
  
  dumxx = DT * dum_ljmain * (old_xx[j] - ret_x(old_xx[j] ,old_xx[l]));
  dumyy = DT * dum_ljmain * (old_yy[j] - ret_y(old_yy[j] ,old_yy[l])); 
  dumzz = DT * dum_ljmain * (old_zz[j] - old_zz[l]);	 
  
  xx[j] += dumxx;
  yy[j] += dumyy;
  zz[j] += dumzz;
  xx[l] -= dumxx;
  yy[l] -= dumyy;
  zz[l] -= dumzz;
}

void interac_supra_others(int j, int l, double old_xx[], double old_yy[], double old_zz[], double old_r[], double old_agek[], double xx[], double yy[], double zz[])
{
  double distlj, LJ1, dum_lj6, dum_ljmain, dumxx, dumyy, dumzz, spring_force;

  distlj = dist3(old_xx[j], old_yy[j], old_zz[j], old_xx[l], old_yy[l], old_zz[l]);
  
  if (distlj < old_r[j] + old_r[l]) {
    LJ1 = (old_r[j] + old_r[l]) / distlj;
    dum_lj6 = pow(LJ1, 6);
    dum_ljmain = 4.0*eps_m * dum_lj6 * (dum_lj6-1.0) / (distlj*distlj);
  }
  else { //attractive 
    if(old_agek[j] > THRESH_SP && old_agek[l] > THRESH_SP) 
      spring_force = K_TOTAL;
    else
      spring_force = K_DESMOSOME;
    dum_ljmain = adhesion(distlj, old_r[j], old_r[l], spring_force);
  }
  
  dumxx = DT * dum_ljmain * (old_xx[j] - ret_x(old_xx[j] ,old_xx[l]));
  dumyy = DT * dum_ljmain * (old_yy[j] - ret_y(old_yy[j] ,old_yy[l])); 
  dumzz = DT * dum_ljmain * (old_zz[j] - old_zz[l]);	 

  xx[j] += dumxx;
  yy[j] += dumyy;
  zz[j] += dumzz;
  xx[l] -= dumxx;
  yy[l] -= dumyy;
  zz[l] -= dumzz;

}
void interac_lj(int j, int l, double old_xx[], double old_yy[], double old_zz[], double old_r[], 
			double xx[], double yy[], double zz[])
{
  double distlj, LJ1, dum_lj6, dum_ljmain, dumxx, dumyy, dumzz;

  distlj = dist3(old_xx[j], old_yy[j], old_zz[j], old_xx[l], old_yy[l], old_zz[l]);
  
  LJ1 = (old_r[j] + old_r[l]) / distlj;
  dum_lj6 = pow(LJ1, 6);
  dum_ljmain = 4.0*eps_m * dum_lj6 * (dum_lj6-1.0) / (distlj*distlj);

  dumxx = DT * dum_ljmain * (old_xx[j] - ret_x(old_xx[j] ,old_xx[l]));
  dumyy = DT * dum_ljmain * (old_yy[j] - ret_y(old_yy[j] ,old_yy[l])); 
  dumzz = DT * dum_ljmain * (old_zz[j] - old_zz[l]);	 

  xx[j] += dumxx;
  yy[j] += dumyy;
  zz[j] += dumzz;
  xx[l] -= dumxx;
  yy[l] -= dumyy;
  zz[l] -= dumzz;

}


void interac_stem_membrane(int j, int l, double old_xx[], double old_yy[], double old_zz[], double old_r[], 
			   int old_state[], double spring_force, double xx[], double yy[], double zz[])
{
  double dum_ljmain, distlj, LJ1, LJ2, dum_lj6, dumxx, dumyy, dumzz;
  
  dum_ljmain = 0.0;
  distlj = dist3(old_xx[j], old_yy[j], old_zz[j], old_xx[l], old_yy[l], old_zz[l]);
  
  if (distlj < old_r[j] + old_r[l]) { // repulsive
    LJ1 = (old_r[j] + old_r[l]) / distlj;
    dum_lj6 = pow(LJ1, 6);
    dum_ljmain = 4.0*eps_m * dum_lj6 * (dum_lj6-1.0) / (distlj*distlj);
  }
  else { // attractive 
    // in case that coupled musumes still stay in the basal layer even after they lost dermal cell
    if (old_state[j] == FIX) {
      LJ2 = distlj / (old_r[j] + old_r[l]);
      dum_ljmain = -(spring_force/distlj) * (LJ2-1.0);
    }
    else  dum_ljmain = adhesion(distlj, old_r[j], old_r[l], spring_force);
  }
  if (fabs(dum_ljmain) > 100.0) {
    printf("error: interaction stem-%d too strong\n", old_state[l]);
    printf("id=%d xx=%f yy=%f zz=%f ljmain=%f\n", j, xx[j], yy[j], zz[j], dum_ljmain);
    exit(1);
  }
  
  dumxx = DT * dum_ljmain * (old_xx[j] - ret_x(old_xx[j] ,old_xx[l]));
  dumyy = DT * dum_ljmain * (old_yy[j] - ret_y(old_yy[j] ,old_yy[l])); 
  dumzz = DT * dum_ljmain * (old_zz[j] - old_zz[l]);	 
  xx[j] += dumxx;
  yy[j] += dumyy;
  zz[j] += dumzz;
  xx[l] -= dumxx;
  yy[l] -= dumyy;
  zz[l] -= dumzz;
}

void interac_stem_stem(int j, int l, double old_xx[], double old_yy[], double old_zz[], double old_r[], 
		       int other_cell[], double xx[], double yy[], double zz[])
{
  double dum_ljmain, distlj, LJ1, dum_lj6, dumxx, dumyy, dumzz;

  if (other_cell[j] == l) dum_ljmain = 0.0;
  else {
    distlj = dist3(old_xx[j], old_yy[j], old_zz[j], old_xx[l], old_yy[l], old_zz[l]);
    
    if (distlj < old_r[j] + old_r[l]) { // repulsive
      LJ1 = (old_r[j] + old_r[l])/distlj;
      dum_lj6 = pow(LJ1, 6);
      dum_ljmain = 4.0*eps_m * dum_lj6 * (dum_lj6-1.0) / (distlj*distlj);
    }
    else  // adhesion 
      dum_ljmain = adhesion(distlj, old_r[j], old_r[l], K_DESMOSOME);
  }
  
  if (fabs(dum_ljmain) > 100.0) {
    printf("error: interaction stem-stem too strong\n");
    printf("id=%d xx=%f yy=%f zz=%f ljmain=%f %d\n", j, xx[j], yy[j], zz[j], dum_ljmain, l);
    exit(1);
  }
  
  dumxx = DT * dum_ljmain * (old_xx[j] - ret_x(old_xx[j] ,old_xx[l]));
  dumyy = DT * dum_ljmain * (old_yy[j] - ret_y(old_yy[j] ,old_yy[l])); 
  dumzz = DT * dum_ljmain * (old_zz[j] - old_zz[l]);	 
  xx[j] += dumxx;
  yy[j] += dumyy;
  zz[j] += dumzz;
  xx[l] -= dumxx;
  yy[l] -= dumyy;
  zz[l] -= dumzz;
}

void setup_memb_indices(int jr[], int jl[], int jll[], int ju[], int jb[], int jbb[])
{
  int j, jj, kk;
  for (j=0; j<NMEMB; j++) {
    jj = j % NMX;
    kk = j / NMX;
    if (jj == 0) jl[j] = j + NMX-1;
    else jl[j] = j-1;
    
    if (jj <= 1) jll[j] = j + NMX-2;
    else jll[j] = j-2;
    
    if (jj == NMX-1) jr[j] = j - (NMX-1);
    else jr[j] = j+1;
    
    if (kk == 0) jb[j] = j + NMX*NMY - NMX;
    else jb[j] = j-NMX;
    
    if (kk <= 1) jbb[j] = j + NMX*NMY - 2*NMX;
    else jbb[j] = j-2*NMX;
    
    if (kk == NMY-1) ju[j] = j - (NMX*NMY - NMX);
    else ju[j] = j + NMX; 
    
    if (jl[j] < 0 || jl[j] >= NMEMB) { printf("error: illegal index\n"); exit(1);}
    if (jll[j] < 0 || jll[j] >= NMEMB) { printf("error: illegal index\n"); exit(1);}
    if (jr[j] < 0 || jr[j] >= NMEMB) { printf("error: illegal index\n"); exit(1);}
    if (jb[j] < 0 || jb[j] >= NMEMB) { printf("error: illegal index\n"); exit(1);}
    if (jbb[j] < 0 || jbb[j] >= NMEMB) { printf("error: illegal index\n"); exit(1);}
    if (ju[j] < 0 || ju[j] >= NMEMB) { printf("error: illegal index\n"); exit(1);}
  }
}

void bend_interac(double old_xx[], double old_yy[], double old_zz[], double ipn[], double ipm[],
		  double dn[], double dm[], double nx[], double ny[], double nz[], double mx[], double my[], double mz[], 
		  int jr[], int jl[], int jll[], int ju[], int jb[], int jbb[], double *distmax, double *distmin)
{
  int j;

  for (j=0; j<NMEMB; j++) {
    dn[j] = dist3(old_xx[jr[j]], old_yy[jr[j]], old_zz[jr[j]], old_xx[j], old_yy[j], old_zz[j]);

    if (dn[j] < 1.0e-5) {
      printf("error: distance dn too short\n");
      exit(1);
    }

    nx[j] = perio_diff(old_xx[jr[j]], old_xx[j], LX) / dn[j]; 
    ny[j] = perio_diff(old_yy[jr[j]], old_yy[j], LY) / dn[j]; 
    nz[j] = (old_zz[jr[j]] - old_zz[j]) / dn[j];

    dm[j] = dist3(old_xx[ju[j]], old_yy[ju[j]], old_zz[ju[j]], old_xx[j], old_yy[j], old_zz[j]);

    if (dm[j] < 1.0e-5) {
      printf("error: distance dm too short\n");
      exit(1);
    }

    mx[j] = perio_diff(old_xx[ju[j]], old_xx[j], LX) / dm[j]; 
    my[j] = perio_diff(old_yy[ju[j]], old_yy[j], LY) / dm[j]; 
    mz[j] = (old_zz[ju[j]] - old_zz[j]) / dm[j];
/*
    if (*distmax < dn[j]) *distmax = dn[j];
    else if (*distmin > dn[j]) *distmin = dn[j];
    if (*distmax < dm[j]) *distmax = dm[j];
    else if (*distmin > dm[j]) *distmin = dm[j];
*/
  }

  for (j=0; j<NMEMB; j++) {
    ipn[j] = nx[jr[j]] * nx[j] + ny[jr[j]] * ny[j] + nz[jr[j]] * nz[j];
    ipm[j] = mx[ju[j]] * mx[j] + my[ju[j]] * my[j] + mz[ju[j]] * mz[j];
  }
}

void interac_wall(double old_zz, double r, double *zz) 
{
  double distlj, LJ1, dum_lj6, dum_ljmain;
  distlj = 2*old_zz;
  double d_inv = 1/distlj;
  LJ1 = 2*r*d_inv;
  dum_lj6 = pow_hexa(LJ1);
  dum_ljmain = 4.0*eps_m * dum_lj6 * (dum_lj6-1.0) *d_inv*d_inv ;
      
  *zz += DT * dum_ljmain * 2.0 * old_zz;	 
}



void interac_memb_memb(int j, int l, double old_xx[], double old_yy[], double old_zz[], double old_r[], 
		       int old_state[], double xx[], double yy[], double zz[])
{
  double distlj, cr_dist, lambda_dist, dum_lj6, dum_ljmain, dumxx, dumyy, dumzz;
  //ra=old_r[j]; rb=old_r[l];
  double rad_sum = old_r[j]+old_r[l];
  double diffx = perio_diff_x(old_xx[j],old_xx[l]);
  double diffy = perio_diff_x(old_yy[j],old_yy[l]);
  double diffz = old_zz[j] - old_zz[l];
  double distljSq = diffx*diffx+diffy*diffy+diffz*diffz;
  //distlj = sqrt(diffx*diffx+diffy*diffy+diffz*diffz);
  //distlj = dist3(old_xx[j], old_yy[j], old_zz[j], old_xx[l], old_yy[l], old_zz[l]);

  cr_dist = rad_sum * P_MEMB;
  double cr_distSq=cr_dist*cr_dist;

  
  if (distljSq < cr_distSq) {
    dum_lj6 = cr_distSq*cr_distSq*cr_distSq/(distljSq*distljSq*distljSq);
    dum_ljmain = 4*eps_m * dum_lj6*(dum_lj6-1.0) / distljSq;
  }
  else if (distljSq < rad_sum*rad_sum){
      //printf("srrtt\n");
      distlj=sqrt(distljSq);
    dum_ljmain = DER_DER_CONST/distlj-DER_DER_CONST/cr_dist;//-(DER_DER_CONST/distlj) * (distlj/cr_dist - 1.0);
  }else {
      lambda_dist = (1.0+P_MEMB)*rad_sum;
        distlj=sqrt(distljSq);
    dum_lj6 = pow_hexa(cr_dist / (lambda_dist - distlj));
    dum_ljmain = -(DER_DER_CONST/ rad_sum)*((1.0-P_MEMB)/P_MEMB)
      -4*eps_m*(dum_lj6*(dum_lj6-1.0)) / ((lambda_dist - distlj)*distlj);
  }
/*
  if (fabs(dum_ljmain) > 100.0) {
    printf("error: interaction der-der too strong\n");
    printf("id=%d state=%d l=%d state=%d xx=%f yy=%f zz=%f, dum_ljmain=%f\n", j, old_state[j], l, old_state[l], xx[j], yy[j], zz[j], dum_ljmain);
    exit(1);
  }
*/
  dumxx = DT * dum_ljmain * diffx;
  dumyy = DT * dum_ljmain *diffy;
  dumzz = DT * dum_ljmain * diffz;
  
  xx[j] += dumxx;
  yy[j] += dumyy;
  zz[j] += dumzz;
  xx[l] -= dumxx;
  yy[l] -= dumyy;
  zz[l] -= dumzz;
}

#define delta_lipid 0.05 // 1.0
double k_lipid(double u , double age)
{
  return 0.25*para_lipid*(1+tanh((ubar-u)/0.01))*(1 + tanh((age-THRESH_SP)/delta_lipid));
}

double k_lipid_release(double u, double age)
{
    return 0.25*para_lipid_rel*(1 + tanh((u - ubar)/0.01))*(1 + tanh((age  -  THRESH_SP)/delta_lipid));
}


double Kx(double age, int state)
{
  if(state == ALIVE || state == DEAD){
    return  1.0; // + ((1.0 - 0.1) / 2.) * (1.0 + tanh((22. - age) / 0.001));
  }
  else {
    return 1.0;
  } 
}

double Ky(double age, int state)
{
  if(state == ALIVE || state  == DEAD){
    return  1.0; // + ((1.0 - 0.1) / 2.) * (1.0 + tanh((22. - age) / 0.001));
  }
  else {
    return 1.0;
  }
}

double fxy(double x, double y ){
  double x1,y1,myux,myuy,myu,height,sigma;
  myu=12.5;
  height=900.0;
  sigma=4.0;
  if(y<=LY/2&&x<=LX/2){
    myux=myu;
    myuy=myu;
    
  }else if(y>=LY/2&&x>=LX/2){
    myux=myu+LX/2;
    myuy=myu+LY/2;
  }else{
    myux=myu;
    myuy=myu;
  }
  
  x1=exp(-pow((x-myux),2)/(2.0*pow(sigma,2)))/(sqrt(2.0*M_PI)*sigma);
  y1=exp(-pow((y-myuy),2)/(2.0*pow(sigma,2)))/(sqrt(2.0*M_PI)*sigma);
  return height*x1*y1;
}
double Dfx(double x,double y,double z){
	double x1,y1,myux,myuy,myu,height,sigma;
	myu=12.5;
	height=900.0;
	sigma=4.0; if(y<=LY/2){
		myux=myu;
	}else{
		myux=myu+LENGTH_X/2;
	}
	return fxy(x,y)*(-x+myux)/pow(sigma,2);
}    
double Dfy(double x,double y,double z){
	double x1,y1,myux,myuy,myu,height,sigma;
	myu=12.5;
	height=900.0;
	sigma=4.0;
        if(y<=LY/2){
		myuy=myu;
	}else{
		myuy=myu+LENGTH_Y/2;
	}
	return fxy(x,y)*(-y+myuy)/pow(sigma,2);

}
double Dfz(double x,double y,double z){
	return 1.0;
}

double f(double norm, double r1, double r2)
{
  return para_k * (0.95*(r1 + r2) - norm) / norm / para_mu;
}

double ageb_const(int state , double u, int tb)
{
  double amp;
  if (tb < MALIGNANT) amp = accel_div;
  else amp = 1.0;

  return amp * eps_kb * (u0 + para_alpha_b * positive(u - u0));


}

double agek_const(int state , double u, int tb)
{
  double amp;
  if (state == DEAD || state == AIR) 
    return eps_kk * u0;
  else if (state == ALIVE) {
    if (tb < MALIGNANT) amp = accel_diff;
    else amp = 1.0;
    
    return amp * eps_ks * (S0 + para_alpha_k * positive(u - u0));
  }
  else {
    printf("error: state must be DEAD (AIR) or ALIVE\n");
    exit(1);
  }
}

double fr(double r, int state)
{
  if(r < R_max)
    return para_er * (R_max - r) * r;
  else return 0;
}

void find_dermis(int j, double xx[], double yy[], double zz[], int state[], int lj[NN][N2], int *p1, int indx[])                                           
{                                                                            
  int k;
  int l;
  double d1 = LENGTH_X, d2 = LENGTH_X, d3 = LENGTH_X;                        
  double distance;                                                           
  *p1 =  - 1;                                                                
  for (k = 0; k < indx[j]; k++) {                                            
    l = lj[j][k];                         
    if(state[l] == MEMB){                                                     
      distance = dist3(xx[l], yy[l], zz[l], xx[j], yy[j], zz[j]);            
      if(distance < d1){                                                     
        d1 = distance;                                                       
        *p1=l;                                                               
      }                                                                      
    }                                                                        
  }                                                                          
}   

void div_direction(int j, double xx[], double yy[], double zz[], double r[], int old_state[], 
		   int indx[], double *nx, double *ny, double *nz, int p1)
/* modified: 
   dermal cell label p1 is calculated and supplied from outside */
{
  double rand_arg, rand_theta, rand_phi, cr1, sr1, cr2, sr2;
  double randnx, randny, randnz, sum;
  double mx, my, mz;
  double lx, ly, lz;

  calc_normal(j, p1, xx, yy, zz, r, &lx, &ly, &lz, &mx, &my, &mz);

  // prepare a random unit vector 
  do {
    rand_theta=M_PI*genrand_res53();
    cr1 = cos(rand_theta);
    sr1 = sin(rand_theta);
    rand_phi=2*M_PI*genrand_res53();
    cr2 = cos(rand_phi);
    sr2 = sin(rand_phi);
    
    randnx = sr1 * cr2;
    randny = sr1 * sr2;
    randnz = cr1;
    
    // take the outer product of the normal l and a random unit vector
    *nx = ly * randnz - lz * randny;
    *ny = lz * randnx - lx * randnz;
    *nz = lx * randny - ly * randnx;
  }  
  while((sum = sqrt((*nx)*(*nx) + (*ny)*(*ny) + (*nz)*(*nz))) < 1.0e-7);

  *nx /= sum;
  *ny /= sum;
  *nz /= sum;

}

void set_newstate(int *oldst, int *st){
  /* if symmetric */
  //  *oldst = *st = FIX;

  /* if assymmetric */
  *oldst = *st = MUSUME;
}

void calc_normal(int i, int p1, double x[],double y[],double z[], double r[], double *nx, double *ny, double *nz, double *mx, double *my, double *mz){
  double nx1, ny1, nz1;
  double bx1, by1, bz1;
  double bx2, by2, bz2;
  double d, distance, cos_theta;

  if (p1==-1) {
    printf("error in function calc_normal: dermis not found\n");
    exit(1);
  }
  
  nx1 = x[i] - ret_x(x[i], x[p1]);
  ny1 = y[i] - ret_y(y[i], y[p1]);
  nz1 = z[i] - z[p1];
  
  *nx  =  nx1/norm(nx1, ny1, nz1);
  *ny  =  ny1/norm(nx1, ny1, nz1);
  *nz  =  nz1/norm(nx1, ny1, nz1);
  
  distance = dist3(x[i],y[i],z[i],x[p1],y[p1],z[p1]);
  
  *mx = x[i] - distance * (*nx);
  *my = y[i] - distance * (*ny);
  *mz = z[i] - distance * (*nz);
}

double adhesion(double distlj, double rj, double rl, double spring_const)
{
  double LJ2;

  LJ2 = distlj / (rj + rl);
  return (LJ2 > LJ_THRESH ? 
	  0.0 : 
	  -(spring_const/distlj) 
	  * ((LJ2-1.0) - pow(LJ2-1.0,3)/((LJ_THRESH-1.0)*(LJ_THRESH-1.0))));
}
