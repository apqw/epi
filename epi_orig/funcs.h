#include <stdio.h>

/*********-- in 21-der.c **************/
void output_last_data(int ncell, double xx[], double yy[], double zz[],  
		      int state[], double r[], double ageb[], double agek[], double u_ave[], 
		      int div_times[], double fat[], double nut[], int touch[], int other_cell[],int tb[], double u[], double v[], 
		      double p[], double w[NN][N2], double a[NX+1][NY+1][NZ+1], double B[NX+1][NY+1][NZ+1]);
void cell_dynamics(int *state, double *ageb, double *agek, 
		   double *xx, double *yy, double *zz, double *r, 
		   int lj[NN][N2], int indx[], 
		   double *u, double *v, double *p, double *fat, double *Vc, 
		   double w[NN][N2], double *u_ave, int *div_times, 
		   int *ncell, int *nvac, int vanished_cell[], int i, int *sw, 
		   int *nborn, int *ndisap, int *touch, int *pair, int *pair2, int *pair_indx, 
		   double L[NN], int other_cell[NN], int ljd[NN][N2],int ljd_indx[], int tb[NN]);
void calc_b(double old_B[NX+1][NY+1][NZ+1], int *old_state, double *old_agek, 
	    int map[NX+1][NY+1][NZ+1], int map2[NX+1][NY+1][NZ+1], double zzmax, 
	    double B[NX+1][NY+1][NZ+1]);
void ca_dynamics(double *u_ave, double B[NX+1][NY+1][NZ+1], 
		 double *u, double *v, double *p, double a[NX+1][NY+1][NZ+1], 
		 double w[NN][N2], int *ixx, int *iyy, int *izz, int gj[NN][N2], int indx[],
		 int *state, double *ageb, double *agek, int map[NX+1][NY+1][NZ+1], 
		 int map2[NX+1][NY+1][NZ+1], double zzmax, 
		 double *r, double *xx, double *yy, double *zz, int n);

void evolution(double u[NN], double v[NN], double p[NN], double w[NN][N2], double a[NX+1][NY+1][NZ+1],
	  double B[NX+1][NY+1][NZ+1], double xx[NN], double yy[NN], double zz[NN], 
	  double r[NN], double ageb[NN], double agek[NN], 
	       int state[NN], int div_times[NN], double fat[NN], double Vc[NN], int touch[NN],
	       int lj[][N2], int indx[], int *ncell_var, double u_ave[NN], int pair[NN], int pair2[NN], int *pair_indx, double L[NN], int other_cell[NN], int ljd[NN][N2], int ljd_indx[], int tb[]); 


/************* func21-2.c ******************/
void initial_u(double u[], double v[], double p[], double a[NX+1][NY+1][NZ+1], double B[NX+1][NY+1][NZ+1], double w[NN][N2]);
void save_start(FILE **fp, char *file_name);
void save_start_a(FILE **fp, char *file_name);
void save_end(FILE **fp);
void check_pair(int ncell, double xx[], double yy[], double zz[], 
		int pair[], int pair2[], int *pair_indx, double L[], int other_cell[]);
void connect_lj(int lj[][N2], int state[], int ncell, double xx[], double yy[], double zz[], double r[], int indx[]);
void input_conc(double u[], double v[], double p[], double w[][N2], double a[NX+1][NY+1][NZ+1], double B[NX+1][NY+1][NZ+1]
		, double ini_u[], double ini_v[], double ini_p[], double ini_w[][N2], 
		double ini_a[NX+1][NY+1][NZ+1], double ini_B[NX+1][NY+1][NZ+1]);
void check_zzmax(double zz[], int state[], int ncell, double *zzmax);
void initialize_vars(FILE *fi, double xx[], double yy[], double zz[], 
		     int lj[NN][N2], int state[], double r[], double ageb[], double agek[], int div_times[], 
		     double fat[], double Vc[], int touch[], double L[],
		     int other_cell[], int tb[], int *ncell_var, int *nder, int *nmemb);
void reopen_file(FILE **f, char *str);

void output_data(int i, double xx[], double yy[], double zz[], int lj[NN][N2],int state[], double r[], 
		 double ageb[], double agek[], double u[], double w[NN][N2], double av[], 
		 int div_times[], double L[], double fat[],double Vc[],int touch[], int other_cell[], int tb[],
                 int indx[], int ncell, double a[NX+1][NY+1][NZ+1], 
		 double B[NX+1][NY+1][NZ+1]);
void initialize_state(int *state, double *ageb, double *agek, double *fat, 
		      double *Vc);

void check_area(int *old_state, double *old_xx, double *old_yy, double *old_zz,
		int area[ANX][ANY][ANZ][N3], int aindx[ANX][ANY][ANZ], int n);
void setup_map(int map[NX+1][NY+1][NZ+1], int map2[NX+1][NY+1][NZ+1], int state[], int ixx[], int iyy[], int izz[], 
	       double xx[], double yy[], double zz[], double r[], int n);
void output_data_on_map(FILE **fp2, FILE **fp3, int map[NX+1][NY+1][NZ+1], 
			int state[], double u[], double v[], double p[], 
			double agek[], double fat[], 
			double a[NX+1][NY+1][NZ+1],double B[NX+1][NY+1][NZ+1], 
			int kitei_c, int i);
void convert_coordinate_to_lattice(double xx[], double yy[], double zz[], 
				   int ixx[], int iyy[], int izz[], int n);
void errorcheck(int ncell);
void check_vacant_seats(int ncell, int *ndisap, int vanished_cell[], int state[]);
void calc_energy(int i, double xx[], double yy[], double zz[], int state[], int div_times[],
		 double u[], int touch[], int tb[], int ncell, int nborn, int ndisap, int NFIX);
void check_localization(int lj[NN][N2], int indx[], int state[], int touch[], int ncell);
void print_data(int time, int param_div_max, int num_n0);
void input_others(double u_ave[], double u[], double v[], double p[], double r[], double xx[], 
		  double yy[], double zz[], double ageb[], double agek[]
		  , int state[], int div_times[], int ncell);
void recycle(int n,double xx[], double yy[], double zz[], int state[], double r[], 
		 double ageb[], double agek[], double av[], 
	     int div_times[], double L[], double fat[],double Vc[],int touch[], 
	     int other_cell[],int tb[],  double u[],double v[],double p[],
	     double w[NN][N2],double a[NX+1][NY+1][NZ+1],double B[NX+1][NY+1][NZ+1]);

void initialize_sc(int state[], double zz[], double zzmax, double agek[], int ncell);

double ret_x(double x1, double x2);
double ret_y(double x1, double x2);
double norm(double x, double y, double z);
double dist(double x1, double y1, double x2, double y2);
double dist3(double x1, double y1, double z1, double x2, double y2, double z2);
double dist3_inv(double x1, double y1, double z1, double x2, double y2, double z2);
inline double pow_hexa(double x){
    double sq = x*x;
    return sq*sq*sq;
}
inline double normSq(double x,double y,double z){
    return x*x+y*y+z*z;
}
inline double dist3Sq(double x1, double y1, double z1, double x2, double y2, double z2){
    return normSq(x1-ret_x(x1,x2),y1-ret_y(y1,y2),z1-z2);
}

int count_line(FILE *);

/*----- in ca_dynamics -----*/
double J(double p);
double fu(double u, double v, double p, double a, double B, int debug_j, int ii);
double fv(double u, double v, double p, double a, double age, int state);
double fp(double u, double v, double p, double a, double age, int state);
double IAG(double age, int state);
double fw(double diff, double w);
double fa(double u, double a);
double I(double u);
double g(double a);

void setup_air_stim(int map[NX+1][NY+1][NZ+1], int lj[][N2], int indx[], int state[], double air_stim[NX+1][NY+1][NZ+1]);

/*----- in cell_dynamics -----*/
double Kx(double age, int state);
double Ky(double age, int state);
double fxy(double x, double y );
double Dfx(double x,double y,double z);
double Dfy(double x,double y,double z);
double Dfz(double x,double y,double z);
double f(double norm, double r1, double r2);
double positive(double x);
double ageb_const(int state , double u, int tb);
double agek_const(int state , double u, int tb);
double fr(double r, int state);
double fB(double age, double B, int state);
double IR(double age, int is_cornified);
double k_lipid(double u , double age);
double k_lipid_release(double u, double age);

void div_direction(int j, double xx[], double yy[], double zz[], double r[], int old_state[], 
		   int indx[], double *nx, double *ny, double *nz, int p1);
void div_direction_psoriasis(int j, double xx[], double yy[], double zz[], double r[], int old_state[], 
		   int lj[][N2], int indx[], double *nx, double *ny, double *nz);
void div_direction_2D(int j, double xx[], double yy[], double zz[], double r[], int old_state[], 
		   int lj[][N2], int indx[], double *nx, double *ny, double *nz);

void set_newstate(int *oldst, int *st);
void check_differentiation_disa(int ncell, int state[], double zz[], int other_cell[], int pair_indx);

int is_within_cell(double x1, double y1, double z1, double x2, double y2, double z2, double r);
inline int is_within_cell_sq(double x1, double y1, double z1, double x2, double y2, double z2, double r){
    return (dist3Sq(x1, y1, z1, x2, y2, z2)<r*r);
}

void calc_radius(double *ra, double *rb, double *dum_dist, double x1, double y1, double z1, double r1, 
		 double x2, double y2, double z2, double r2);
void find_dermis(int j, double xx[], double yy[], double zz[],int state[], int lj[NN][N2], int *p1, int indx[]);
void calc_normal(int i, int p1, double x[],double y[],double z[], double r[], double *nx, double *ny, double *nz, double *mx, double *my, double *mz);

void check_lj_connection(double xx[], double yy[], double zz[], int lj[][N2],int indx[], int ncell);
double DU_dermis_spring(double dist, double r);
void divide_cell(int j, double xx[], double yy[], double zz[], double r[], 
		 int old_state[], int state[], int div_times[], double u[], double u_ave[], 
		 double v[], double p[], double old_ageb[], double old_agek[], 
		 double ageb[], double agek[], double old_fat[], double fat[], double old_Vc[], double Vc[], 
		 int *ncell, int *nvac, int *nborn,
		 int vanished_cell[], int divtimes_inherited, 
		 int i, int lj[][N2], int indx[], int *pair, int *pair2, 
		 int *pair_indx, double L[NN], int other_cell[], int tb[], int p1);
void periodic_cond(double xx[], double yy[], int ncell);
void initialize_old_values(double *old_xx, double *old_yy, double *old_zz, int *old_state, 
			   double *old_ageb, double *old_agek, double *old_r,  
			   double *old_fat, double *old_Vc, 
			     double *xx, double *yy, double *zz, int *state, 
			   double *ageb, double *agek, double *r, 
			   double *fat, double *Vc, int ncell);
void moveover(int *pair_indx, int pair[], int pair2[], int nunpaired);
double bend_force(double n[], double m[], double ipn[], double ipm[], double dn[], double dm[], 
	     int j, int jr, int jl, int jll, int ju, int jb, int jbb);
double bend_force_sqr(double n[], double m[], double ipn[], double ipm[], double dn[], double dm[], 
	     int j, int jr, int jl, int jll, int ju, int jb, int jbb);
double perio_diff(double x1, double x2, double lb);
inline double perio_diff_x(double x1, double x2)
{
double diff=x1-x2;
  if (diff*diff < 0.25*LX*LX) return diff;
  else {
    if (x1 > x2) {
      return diff-LX;
    }
    else {
      return diff+LX;
    }
  }
}

inline double perio_diff_y(double x1, double x2)
{
double diff=x1-x2;
  if (diff*diff < 0.25*LY*LY) return diff;
  else {
    if (x1 > x2) {
      return diff-LY;
    }
    else {
      return diff+LY;
    }
  }
}
void set_spring_force(int j, int old_state[], int other_cell[], int div_times[], double *spring_force);
void interac_lj_repulsive(int j, int l, double old_xx[], double old_yy[], double old_zz[], double old_r[], 
			  double xx[], double yy[], double zz[]);
void interac_stem_membrane(int j, int l, double old_xx[], double old_yy[], double old_zz[], double old_r[], 
			   int old_state[], double spring_force, double xx[], double yy[], double zz[]);
void interac_stem_stem(int j, int l, double old_xx[], double old_yy[], double old_zz[], double old_r[], 
		       int other_cell[], double xx[], double yy[], double zz[]);
void setup_memb_indices(int jr[], int jl[], int jll[], int ju[], int jb[], int jbb[]);
void bend_interac(double old_xx[], double old_yy[], double old_zz[], double ipn[], double ipm[],
		  double dn[], double dm[], double nx[], double ny[], double nz[], double mx[], double my[], double mz[], 
		  int jr[], int jl[], int jll[], int ju[], int jb[], int jbb[], double *distmax, double *distmin);
void interac_wall(double old_zz, double r, double *zz) ;
void interac_memb_memb(int j, int l, double old_xx[], double old_yy[], double old_zz[], double old_r[], 
		       int old_state[], double xx[], double yy[], double zz[]);
void interac_der_der(int j, int l, double old_xx[], double old_yy[], double old_zz[], double old_r[], 
		     double xx[], double yy[], double zz[]);
void interac_lj(int j, int l, double old_xx[], double old_yy[], double old_zz[], double old_r[], 
			double xx[], double yy[], double zz[]);
double adhesion(double distlj, double rj, double rl, double spring_const);
void interac_supra_others(int j, int l, double old_xx[], double old_yy[], double old_zz[], double old_r[], double old_agek[], double xx[], double yy[], double zz[]);
