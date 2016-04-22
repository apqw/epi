double debug_double;
int debug_int;

int ca_N = 1000; // 1000 // #steps for Ca averaging 
double para_thgra = 0.2; //done
double para_thpri = 1.0; //done
double delta_th = 1.0; //done
double FMAX=2.0;
double FMIN=0.0;
double a0c = 3.;
double lez1 = (2. * Rad);


double lez = 5 * Rad;
double lbx = LENGTH_X * 0.5;
//double ab = 0.05; //0.02; //0.02;
double p0c = 0.;
double ue = 1.5;//1.5;

double Da = 1.;
double du = 0.01;
double dp = 0.1;
double para_kpa = 4.; //done
double para_wd = 0.1;//0.2; //done
double delta_cos = 0.1;
double para_dpri = 0.;

//double eps_kb = 0.10 ;
double eps_kb = 0.1;
double eps_ks = 0.05;      //細胞の状態に応じて状態変数の変動率を変える？
double eps_kk = 0.05;

double B0 = 0.0;

double DB = 0.0009;    //Bの計算で用いる変化率
double delta_IR = 3.;
double para_kb = 0.025;

double para_kbc = 0.4; //done
double para_Hb = 0.01; //do

double para_Cout = 1.0;

/*-------------------*/
double para_Kfat = 0.05;
double calory = 2.0;

double myu = 0.05;

double GF=0.003;

double dx = LX/(double)NX;
double dy = LY/(double)NY;
double dz = LZ/(double)NZ;

int irx = (int) (Rad*NX/LENGTH_X);
int iry = (int) (Rad*NY/LENGTH_Y);
int irz = (int) (Rad*NZ/LENGTH_Z);

double xc, yc, zc;

double LJ_THRESH = 1.2;

int NDER, NMEMB, NTRI;

double AIR_STIM = 0.1;
double accel_div = 1.0; // normal: 1.0
double accel_diff = 1.0; // normal: 1.0
int MALIGNANT = 0; // number of malignant stem cells
double DD_FAC;
double DER_DER_CONST = 0.2;
double PRESS_CONST;
double DD_LJ_FAC;
double KBEND = 0.5;
double distmax=0.0, distmin=0.0;
int SYSTEM = WHOLE;
int KEY_DERMAL_CHANGE = 1;
