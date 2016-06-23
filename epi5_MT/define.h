#ifndef DEFINE_H
#define DEFINE_H
//#define NDEBUG
#include <string>
#include <memory>
#define cdefd static constexpr double
#define cdefi static constexpr int
#define cdefui static constexpr unsigned int
#define cdefs static constexpr size_t

#define defd static double
#define rdefd static double&
#define defi static int
#define rdefi static int&
#define defui static unsigned int
#define rdefui static unsigned int&
#define defs static size_t
#define rdefs static size_t&

#ifdef _WIN32
#define RESTRICT __restrict
#else
#define RESTRICT __restrict__
#endif

class Cell;
using CellPtr = Cell*;//std::shared_ptr<Cell>;

using cmask_ty = uint_fast8_t;
using CELL_STATE_internal = uint_fast8_t;

/*
	DO NOT CHANGE 
*/
enum CELL_STATE:CELL_STATE_internal {
    ALIVE = 0,
    DEAD = 1,
    DISA = 2,
    UNUSED = 3, //
    FIX = 4,
    BLANK = 5,
    DER = 6,
    MUSUME = 7,
    AIR = 8,
    MEMB = 9
};

enum CELL_STATE_MASK :uint_fast16_t {
	ALIVE_M = 1u<<ALIVE,
	DEAD_M = 1u << DEAD,
	DISA_M = 1u<<DISA,
	UNUSED_M = 1u<<UNUSED, //
	FIX_M = 1u<<FIX,
	BLANK_M = 1u<<BLANK,
	DER_M = 1u<<DER,
	MUSUME_M = 1u<<MUSUME,
	AIR_M = 1u<<AIR,
	MEMB_M = 1u<<MEMB
};

enum BoundaryType {
	DIRICHLET, NEUMANN, PERIODIC
};

namespace cont {
	
	cdefs STATE_NUM = 10;
	cdefs MAX_CONNECT_CELL_NUM = 400;
	cdefs MAX_CELL_NUM = 30000;
	cdefd DT_Cell = 0.01;
	cdefd DT_Ca = 0.01;//0.02;
	cdefi MALIG_NUM = 0;
										 /*
										 LX,LY,LZ:ŒvŽZ—Ìˆæ‚ÌƒTƒCƒY
										 */
	cdefd LX = 50.0;
	cdefd LY = 50.0;
	cdefd LZ = 100.0;

	/*
	NX,NY,NZ:ŒvŽZ—Ìˆæ‚Ì•ªŠ„”
	*/
	cdefs NX = 100;
	cdefs NY = 100;
	cdefs NZ = 200;

	/*
	dx,dy,dz:ŒvŽZ—Ìˆæ‚Ì•ªŠ„•(L**,N**‚©‚çŒvŽZ‚³‚ê‚é)
	*/
	cdefd dx = LX / NX; cdefd inv_dx = NX / LX;
	cdefd dy = LY / NY; cdefd inv_dy = NY / LY;
	cdefd dz = LZ / NZ; cdefd inv_dz = NZ / LZ;

    cdefd COMPRESS_FACTOR = 6;//ok
	
    cdefd THRESH_SP = 1.0; //ok 3.0->1.0
	
	cdefd LJ_THRESH = 1.2;//ok
	
	cdefd R_max = 1.4;//ok
	cdefd R_der = 1.4;//ok
	cdefd R_memb = 1.0;//ok
	
	cdefd THRESH_DEAD = 22.0;//ok

	cdefui NUM_ITR = 4 * ((int)1e6); //twice

    

	cdefd Ca_avg_time = 10.0;
	

	cdefd ca2p_init = 0.122;

	cdefd ex_inert_init = 0.97;

	cdefd IP3_init = 0;

	cdefd gj_init = 0.99;

	cdefd ATP_init = 0;

	cdefd ext_stim_init = 0;

	

	cdefd FAC_MAP = 2.0;//ok
	
	cdefi SW_THRESH = 20;//ok

	cdefd T_TURNOVER = 6000.0;//
	cdefi NUM_SC_INIT = 1;//ok ha/???->19->1
	cdefd stoch_corr_coef = 1;
    
    
    
}

cdefui SYSTEM = 0;
cdefui WHOLE = 0;
cdefui BASAL = 1;
cdefui CUT=1000;
static constexpr auto OUTPUTDIR="output";

static constexpr auto last_data_uvp_name="last_data_uvp";
static constexpr auto last_data_cell_name="last_data_cell";
static constexpr auto last_data_w_name="last_data_w_alt";
static constexpr auto last_data_a_name="last_data_a_alt";
static constexpr auto last_data_B_name="last_data_B_alt";

static constexpr bool STOCHASTIC = true;

static constexpr bool FORCE_CORNIF = true;

template<typename T, size_t X, size_t Y, size_t Z>
class Field;
template<typename T>
using FArr3D = Field<T, cont::NX + 1, cont::NY + 1, cont::NZ + 1>;

template<typename T>
using RawArr3D = T[cont::NX + 1][cont::NY + 1][cont::NZ + 1];
#endif // DEFINE_H
