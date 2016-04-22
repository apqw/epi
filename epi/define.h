#pragma once
#include <immintrin.h>
#include <vector>
#include <cassert>
#define _mm256_full_hadd_ps(v0, v1) \
        (_mm256_hadd_ps(_mm256_permute2f128_ps(v0, v1, 0x20), \
                       _mm256_permute2f128_ps(v0, v1, 0x31)))
#define SET8(val) static const __m256 val##_8
#define SET8s(val) static const __m256 val##_8s
#define SET4d(val) static const __m256d val##_4d

#define SET8_d(cl,val) const __m256 cl::val##_8=_mm256_set1_ps((float)(val))
#define SET8s_d(cl,val) const __m256 cl::val##_8s=_mm256_set1_ps((float)(val))
#define SET4d_d(cl,val) const __m256d cl::val##_4d=_mm256_set1_pd((double)(val))

//Check if the 4-paired-cell contains the cell which has this state 
#define CONTAIN_STATE(_cell,state) (_mm256_testz_pd(_cell.state.smask[(state)], _cell.state.smask[(state)]) == 0)

#define VEC_X_RANGE_VALIDATION(var) do{ assert(var >=0);\
assert(var % 4 == 0);\
}while(0)


typedef double real;	//高速化の際にfloatにした方が良い場合があるかもしれないので実数は一旦realとして定義しておく
						//必要になったらfloatに変える
						//SSEとか使う時には注意
						//floatがdoubleより遅い場合はコンパイル時のオプションをチェック
typedef std::vector<std::vector<std::vector<__m256d>>> _3DScalar4d;
#define DEFC static constexpr //打つの面倒

#define DEFC_VEC(name,val) DEFC real name=(val);SET8s(name);SET4d(name)
#define VEC_INIT(cl,name) SET8s_d(cl,name);SET4d_d(cl,name)

class V3D {

};

//__m256 == float*8
//__m256d == double*4
class VSet4d :V3D {
public:
	__m256d x;
	__m256d y;
	__m256d z;
	//need alignment
	//align出来ないなら_mm256_loadu_ps
	VSet4d() {
		x = _mm256_setzero_pd();
		y = _mm256_setzero_pd();
		z = _mm256_setzero_pd();
	}
	VSet4d(const double* &x_arr, const double* &y_arr, const double* &z_arr) {
		x = _mm256_load_pd(x_arr);
		y = _mm256_load_pd(y_arr);
		z = _mm256_load_pd(z_arr);
	}
	//VSet8& operator+(const VSet8 &a); //avoid creating instances
	void operator+=(const VSet4d &a);
	//VSet8& operator*(const VSet8 &a);
	//VSet8& operator*=(const VSet8 &a);
};

class BoundaryInfo {
public:
	//境界条件enum
	enum BCond {
		NONE, DIRICHLET, NEUMANN, PERIODIC
	};
	//範囲
	int min, max;
	//任意のデータ(へのポインタ)(例えば周期境界条件での反対側の値を与えたりするために使う)
	void* data;

	//境界条件
	BCond min_bcond, max_bcond;

	BoundaryInfo() :min(), max(), min_bcond(), max_bcond() {};
	BoundaryInfo(int _min, int _max, BCond _minbc, BCond _maxbc) :min(_min), max(_max), min_bcond(_minbc), max_bcond(_maxbc) {};
	BoundaryInfo(int _min, int _max) :min(_min), max(_max), min_bcond(NONE), max_bcond(NONE) {};
};

class AREA {

};

class CUBE :AREA {
public:
	BoundaryInfo x, y, z;
};




namespace EPI {

	enum STATE {
		ALIVE = 0,
		DEAD = 1,
		DISA = 2,
		FIX = 4,
		BLANK = 5,
		DER = 6,
		MUSUME = 7,
		AIR = 8,
		MEMB = 9
	};
	typedef STATE CELL_STATE;
	class PState4d {
	public:
		static __m256d all1;
		static __m256d all0;
		__m256d smask[10];
		template<typename T, typename... U>
		__m256d getORMask(T first, U... rest) const; //test
		__m256d getORMask() const;

		template<typename T, typename... U>
		__m256d getANDMask(T first, U... rest) const;
		__m256d getANDMask() const;



	};

	//3次元のベクトルをいくつか(レジスタのサイズに応じて)セットで用意する。

	/*
		定数用クラスはこれを継承しておく
	*/
	class C {
	public:
		DEFC_VEC(zero, 0.0);
		DEFC_VEC(half, 0.5);
		DEFC_VEC(one, 1.0);
		DEFC_VEC(neg_zero, -0.0);
		// http://www.math.utah.edu/~beebe/software/ieee/tanh.pdf
		DEFC_VEC(tanh_s_th, 1.29047814397589243466E-08);
		DEFC_VEC(tanh_l_th, 19.06154746539849600897);
		DEFC_VEC(tanh_m_th, 0.54930614433405404570);
		DEFC_VEC(tanh_dbl_p0, -(0.16134119023996228053E+04));
		DEFC_VEC(tanh_dbl_p1, -(0.99225929672236083313E+02));
		DEFC_VEC(tanh_dbl_p2, -0.96437492777225469787);
		DEFC_VEC(tanh_dbl_q0, 0.48402357071988688686E+04);
		DEFC_VEC(tanh_dbl_q1, 0.22337720718962312926E+04);
		DEFC_VEC(tanh_dbl_q2, 0.11274474380534949335E+03);
		DEFC_VEC(tanh_dbl_q3, 1.00000000000000000000);
		static constexpr int max_cell_num = 3000;
		DEFC_VEC(dt_cell, 0.01);
		DEFC_VEC(dt_ca, 0.02);

		/*
		Lx,Ly,Lz:計算領域の幅
		*/
		DEFC_VEC(Lx, 50.0);
		DEFC_VEC(Ly, 50);
		DEFC_VEC(Lz, 100.0);

		/*
		NX,NY,NZ:分割数
		*/
		static constexpr int NX = 100;
		static constexpr int NY = 100;
		static constexpr int NZ = 200;
		SET8s(NX); SET4d(NX);
		SET8s(NY); SET4d(NY);
		SET8s(NZ); SET4d(NZ);
		/*
		dx,dy,dz:グリッドの幅
		*/
		DEFC_VEC(dx, Lx / (double)NX);
		DEFC_VEC(dy, Ly / (double)NY);
		DEFC_VEC(dz, Lz / (double)NZ);

		DEFC_VEC(dxSq, dx*dx);
		DEFC_VEC(dySq, dy*dy);
		DEFC_VEC(dzSq, dz*dz);

		//実行時に指定するのでconstではない
		static int NMEMB;
		static int NDER;
		static int CELL_START_IDX;
		static int x_prev_arr[NX];
		static int x_next_arr[NX];
		static int y_prev_arr[NY];
		static int y_next_arr[NY];
		static int z_prev_arr[NZ];
		static int z_next_arr[NZ];
	};

	//int current_cell_num;
	class CellSet4d {
	public:

		static constexpr int max_reaction_cell_num = 400;
		__m256d valid_mask;
		VSet4d pos;
		__m256d P; //?
		__m256d c; //?
		__m256d ageb;//?
		__m256d agek;//?
		__m256d diff_c; //dc_dt
		__m256d h;
		PState4d state;
		std::vector<bool> react_flag_4d; //flagged block num must be eq or less than max_reaction_cell_num*4
		std::vector<__m256d> react_mask;
		std::vector<__m256d> w;//?
		CellSet4d() {
			//huge memory usage?
			react_mask.resize(C::max_cell_num / 4);
			w.resize(C::max_cell_num / 4);
			react_flag_4d.resize(C::max_cell_num / 4);
			/*
			react_flag_4dで4つとも相互作用しないCellSet4dを除外してからAVXで計算する。
			最悪でも計算量は通常と同じ、最高で1/4以下
			*/
		};
		void get_lattice(VSet4d& out) const;
		template<typename T, typename... U>
		bool hasState(T s, U... rest) const;
		bool hasState() const;
	};
	class Field_Data {
	public:
		static int current_cell_num;
		static std::vector<CellSet4d> cells;
		static std::vector<CellSet4d> next_cells;
		static _3DScalar4d Ca2P_value;
		static int malignant_stemcell_num;
	};

	/*
		param.hに定義されているものを整理
		extern ~の値はmain.hにある？
	*/
	/*
		Ca2+の所の定数定義
	*/
	class Ca2P :C {
	public:
		/*
			dA,dP,dc,dB:拡散係数
		*/
		DEFC_VEC(dA, 1.0);
		DEFC_VEC(dP, 0.1);
		DEFC_VEC(dc, 0.01); //du
		DEFC_VEC(dB, 0.0009);

		/*
			Kaa,Kpp,Kbb:?
		*/
		DEFC_VEC(Kaa, 0.5);
		DEFC_VEC(Kpp, 0.3);
		DEFC_VEC(Kbb, 0.025);

		/*
			Kac,Kf,Kmu,K1,Kg,Kbc:?
		*/
		DEFC_VEC(Kac, 0.002);//originally STIM11
		DEFC_VEC(Kf, 8.1);
		DEFC_VEC(Kmu, 0.05);
		DEFC_VEC(K1, 0.7);
		DEFC_VEC(Kg, 0.1);
		DEFC_VEC(Kbc, 0.4);

		/*
			mu0,mu1,alpha0,gamma,beta,Hb,H0:?
		*/
		DEFC_VEC(mu0, 0.567);
		DEFC_VEC(mu1, 0.1);
		DEFC_VEC(alpha0, 0.11); //para_b


		//?
		DEFC_VEC(sub1Alpha0, 1.0 - alpha0); //para_bb
		//DEFC real sub1Alpha0 = 1.0 - alpha0; //originally para_bb
		DEFC_VEC(gamma, 2.0);

		DEFC_VEC(beta_zero, 0.02);
		DEFC_VEC(CA_OUT, 1.0);

		DEFC_VEC(beta, beta_zero*CA_OUT);
		DEFC_VEC(Hb, 0.01);
		DEFC_VEC(H0, 0.5);
		DEFC_VEC(Cout, 1.0);

		/*
			K2,wd0,epsw0:?
		*/
		DEFC_VEC(K2, 0.7);
		DEFC_VEC(wd0, 0.1);
		DEFC_VEC(eps_w0, 0.1);//unmagic-numbered
		//epsw0はFwの計算に使う

		/*
			THRESH_DEAD	:dead threshold?
			DUR_DEAD	:dead duration?
			DUR_ALIVE	:alive duration?
			THRESH_SP	:?
		*/
		DEFC_VEC(THRESH_DEAD, 22.0);
		DEFC_VEC(DUR_DEAD, 2.0);
		DEFC_VEC(DUR_ALIVE, 0.5);
		DEFC_VEC(THRESH_SP, 3.0);//age of stratum spinosum ( a.k.a. T_pri)

		/*
			Sk1,Sk2,Ss,tau_g,tau_s,delta_tau,delta_I,delta_k,kg,ks:?
		*/
		DEFC_VEC(Sk1, THRESH_DEAD - DUR_ALIVE); //replace formulas like this by Sk1
		DEFC_VEC(Sk2, THRESH_DEAD + DUR_DEAD);
		DEFC_VEC(Ss, THRESH_SP);
		DEFC_VEC(tau_g, 0.2);
		DEFC_VEC(tau_s, 1.0);
		DEFC_VEC(delta_tau, 1.0);
		DEFC_VEC(delta_I, 1.5);
		DEFC_VEC(delta_k, 1.0);
		DEFC_VEC(kg, 4.0); //para_Kgra=para_kpa
		DEFC_VEC(ks, 6.0);

		//R?
		DEFC real _r = 1;
		DEFC_VEC(_rSq, _r*_r);

		DEFC_VEC(iage_kitei, 0.0);

		//functions

		__m256d G(const VSet4d&, const VSet4d&, const __m256d&);
		__m256d Fc(const CellSet4d&, const __m256d&);
		__m256d FP(const CellSet4d&, const __m256d&);
		__m256d Fh(const CellSet4d&);
		__m256d Fw(const __m256d&, const __m256d&, const __m256d&);
		//__m256d FB(const __m256d&,const __m256d&);
		__m256d tau_h(const CellSet4d&);
		__m256d In(const CellSet4d&);
		__m256d Kpa(const CellSet4d&);

		void refresh_Ca(const CUBE& calc_area,
			const _3DScalar4d& currentCa,
			const std::vector<CellSet4d>& all_cells,
			_3DScalar4d& nextCa);

		void refresh_P_i(int calc_index_min, int calc_index_max,
			const _3DScalar4d& currentCa,
			const std::vector<CellSet4d>& all_current_cells,
			std::vector<CellSet4d>& refreshed_cells);

		void refresh_c_i(int calc_index_min, int calc_index_max,
			const _3DScalar4d& currentB,
			const std::vector<CellSet4d>& all_current_cells,
			std::vector<__m256d>& diff_c_out,
			std::vector<CellSet4d>& refreshed_cells);

		void refresh_h_i(int calc_index_min, int calc_index_max,
			const std::vector<CellSet4d>& all_current_cells,
			std::vector<CellSet4d>& refreshed_cells);

		void refresh_w_i_j(int calc_index_min, int calc_index_max,
			const std::vector<CellSet4d>& all_current_cells,
			std::vector<CellSet4d>& refreshed_cells);
	};

	class CellDiv :C {
		DEFC_VEC(eps_kb, 0.1);
		DEFC_VEC(u0, 0.122);
		DEFC_VEC(alpha_b, 5.0);

		DEFC_VEC(agki_max, 6.0);
		DEFC_VEC(fac, 1.0);
		DEFC_VEC(agki_max_fix, fac*agki_max);
		DEFC_VEC(stoch_div_time_ratio, 0.25);
		
		DEFC_VEC(S_star, agki_max*(1.0 - stoch_div_time_ratio));
		DEFC_VEC(T_musume, agki_max / (eps_kb*u0));
		DEFC_VEC(T_stem_fix, agki_max_fix / (eps_kb*u0));
		DEFC_VEC(T_0_musume, T_musume*(1.0 - stoch_div_time_ratio));
		DEFC_VEC(T_0_stem_fix, T_stem_fix*(1.0 - stoch_div_time_ratio));

		static constexpr bool STOCHASTIC = true;
		void divide_cell
	};

}
__m256d _tanh_poly(const __m256d&);
__m256d tanh_avx(const __m256d&);
__m256d tanh_alt(const __m256d&);
__m256d m256dintmod(const __m256d&, const __m256d&);
__m256d calc_avg8(const VSet4d& lattice_4d, const _3DScalar4d& _3DVal_4d); //no boundary condition
void init();
