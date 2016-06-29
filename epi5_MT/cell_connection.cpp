#include "cell_connection.h"
#include "define.h"
#include "utils.h"
#include "cell.h"
#include "cell_conn_value.h"
#include <cinttypes>

/**
 *  @file 細胞の接続に関する定義
 */

using cint = int_fast16_t;
#define CINT_F SCNiFAST16

/** グリッドサイズ*/
static constexpr double AREA_GRID_ORIGINAL	= 2.0;//ok

/** グリッドサイズ (?)*/
static constexpr double AREA_GRID			= AREA_GRID_ORIGINAL + 1e-7;//ok

/** X方向のグリッド数*/
static constexpr cint	ANX					= (cint)((double)cont::LX / AREA_GRID_ORIGINAL + 0.5);//ok
/** Y方向のグリッド数*/
static constexpr cint	ANY					= (cint)((double)cont::LY / AREA_GRID_ORIGINAL + 0.5);//ok
/** Z方向のグリッド数*/
static constexpr cint	ANZ					= (cint)((double)cont::LZ / AREA_GRID_ORIGINAL);//ok
/** グリッド1つ当たりの細胞格納数上限 */
static constexpr cint	N3					= 200; //max grid cell num //ok
/** 1細胞の接続最大数 */
static constexpr cint	N2					= 400; //max conn num //ok

/**
 *  グリッドの初期化
 *  @param [in] cman CellManager
 *  @param [out] aindx 各グリッドに格納されている細胞の個数
 *  @param [out] area 各グリッドに格納されている細胞への参照
 */
void grid_init(CellManager& cman, std::atomic<cint>(&aindx)[ANX][ANY][ANZ] , Cell* (&area)[ANX][ANY][ANZ][N3]) { //DO NOT RESTRICT ptrs in area
	using namespace cont;
	std::memset(aindx, 0, sizeof(std::atomic<cint>)*ANX*ANY*ANZ);
    cman.all_foreach_parallel_native([&](Cell*const c) {//cannot restrict
		//auto&c = cman[i];
		cint aix, aiy, aiz;
		aix = (cint)((0.5*LX - p_diff_x(0.5*LX, c->x())) / AREA_GRID);
		aiy = (cint)((0.5*LY - p_diff_y(0.5*LY, c->y())) / AREA_GRID);
		aiz = (cint)((min0(c->z())) / AREA_GRID);

        if ((aix >= (cint)ANX || aiy >= (cint)ANY || aiz >= (cint)ANZ || aix < 0 || aiy < 0 || aiz < 0)) {
			printf("err\n");
			printf("cx:%lf cy:%lf cz:%lf\n", c->x(), c->y(), c->z());
            printf("aix:%" CINT_F " aiy:%" CINT_F " aiz:%" CINT_F "\n", aix, aiy, aiz);
			assert(false);
		}

		area[aix][aiy][aiz][aindx[aix][aiy][aiz]++] = c;
		/*
		このassertは消してよい
		*/
		assert(aindx[aix][aiy][aiz] < (cint)N3);
        c->connected_cell.force_set_count(c->state == MEMB ? MEMB_ADJ_CONN_NUM : 0);
	});
}

/**
 *  全細胞間の接続処理
 *  @param [in] cman CellManager
 *  @param [in] aindx 各グリッドに格納されている細胞の個数
 *  @param [in] area 各グリッドに格納されている細胞への参照
 */
void connect_proc(CellManager& cman, const std::atomic<cint>(&aindx)[ANX][ANY][ANZ], Cell*const (&area)[ANX][ANY][ANZ][N3]) {
	
#define LAT_ARR(axis)\
struct lat_arr_##axis {\
	cint idx[AN##axis * 3];\
	lat_arr_##axis() {\
        for (size_t i = 0; i < AN##axis * 3; i++) {\
			idx[i] = i%AN##axis;\
		}\
	}\
	cint operator[](int i) {\
		return idx[i];\
	}\
}
	LAT_ARR(X); LAT_ARR(Y); LAT_ARR(Z);

	static lat_arr_X xidx;
	static lat_arr_Y yidx;
	static lat_arr_Z zidx;
	

	using namespace cont;
    constexpr cint ii = 2;
    constexpr cint iix = ii+1;
	cman.non_memb_foreach_parallel_native([&](Cell*const c) { //cannot restrict
		
		const cint anx = (cint)(c->x() / AREA_GRID);
		const cint any = (cint)(c->y() / AREA_GRID);
		const cint anz = (cint)(c->z() / AREA_GRID);

        assert(!(anx >= (cint)ANX || any >= (cint)ANY || anz >= (cint)ANZ || anx < 0 || any < 0 || anz < 0));
		const size_t my_index = c->get_index();
        const cint xend = anx + iix + (cint)ANX;
		const cint yend = any + ii + (cint)ANY;
		const cint zend = anz + ii + (cint)ANZ;
        for (cint j = anx - iix+(cint)ANX; j <= xend; j++) {
			const cint aix = xidx[j];
            for (cint k = any - ii + (cint)ANY; k <= yend; k++) {
				const cint aiy = yidx[k];
                for (cint l = anz - ii + (cint)ANZ; l <= zend; l++) {
					const cint aiz = zidx[l];
					const cint sz = aindx[aix][aiy][aiz];

					for (cint m = 0; m < sz; ++m) {
						Cell*const o = area[aix][aiy][aiz][m];
						if (my_index <= o->get_index())continue;
						const double diffx = p_diff_x(c->x(), o->x());
						const double diffy = p_diff_y(c->y(), o->y());
						const double diffz = c->z() - o->z();
						const double rad_sum = c->radius + o->radius;
						if (DIST_SQ(diffx, diffy, diffz) <= LJ_THRESH*LJ_THRESH*rad_sum*rad_sum) {
							c->connected_cell.push_back(o);
							o->connected_cell.push_back(c);
							assert(c->connected_cell.size() < N2);
						}
					}
				}
			}
		}
	});
}
/**
 *  GJの値の更新
 */
void gj_refresh(CellManager& cman) {

	cman.all_foreach_parallel_native([](Cell*const RESTRICT c) {

        for(auto&& it = c->gj.begin(),itend=c->gj.end();it!=itend;++it){
            it->second.uncheck();
        }

        c->connected_cell.foreach([&c](Cell*const RESTRICT cptr) {
            c->gj[cptr].check();
        });

        for(auto&& it = c->gj.begin(),itend=c->gj.end();it!=itend;++it){
            if(!it->second.is_checked()){
                it->second.reset();
            }
        }

	});
}

/**
 *  dermisの検索
 *  @attention 接続が済んでいる必要がある。
 */
const Cell* find_dermis(const Cell*const RESTRICT c) {
	double d1Sq = cont::LX*cont::LX;
	double distanceSq = 0;
	const Cell* dermis = nullptr;
    c->connected_cell.foreach([&](const Cell*const RESTRICT cptr) {
		if (cptr->state == MEMB) {
			distanceSq = p_cell_dist_sq(c, cptr);
			if (distanceSq < d1Sq) {//2乗になってなかった
				d1Sq = distanceSq;
				dermis = cptr;
			}
		}
	});
	return dermis;
}

/**
 *  dermisが存在すべき全細胞に対して、dermisへの参照をセット。
 *  (FIXとMUSUMEに対してdermisを検索する。)
 */
inline void set_dermis(CellManager& cman) {
    cman.other_foreach_parallel_native([](Cell*const RESTRICT c) {
		
		if (c->state == FIX || c->state == MUSUME) {
			c->set_dermis(find_dermis(c));
		}
	});
}


/**
 *  細胞の接続に関する処理を行う。
 */
void connect_cell(CellManager & cman)
{
	using namespace cont;
	static std::atomic<cint> aindx[ANX][ANY][ANZ] = {};
	static Cell* area[ANX][ANY][ANZ][N3] = { nullptr };
	

	grid_init(cman,aindx,area);
	connect_proc(cman, aindx, area);
	gj_refresh(cman);
	set_dermis(cman);
	
}
