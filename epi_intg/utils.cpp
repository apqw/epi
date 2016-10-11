#include "utils.h"
#include "global.h"
#include <vector>
std::string operator"" _s(const char* str,std::size_t sz) {
    return std::string(str);
}

real p_diff_x(const real x1, const real x2)
{
    //using namespace cont;
    const real diff = x1 - x2;
    if (diff > 0.5*pm->LX)return diff - pm->LX;
    if (diff <= -0.5*pm->LX)return diff + pm->LX;
    return diff;
}

real p_diff_y(const real y1, const real y2)
{
    //using namespace cont;
    const real diff = y1 - y2;
    if (diff > 0.5*pm->LY)return diff - pm->LY;
    if (diff <= -0.5*pm->LY)return diff + pm->LY;
    return diff;
}

real p_dist_sq(const real x1, const real y1, const real z1, const real x2, const real y2, const real z2)
{
    const real diffx = p_diff_x(x1, x2);
    const real diffy = p_diff_y(y1, y2);
    const real diffz = z1 - z2;
    return diffx*diffx + diffy*diffy + diffz*diffz;
}

real min0(const real a) {
    return a > 0 ? a : 0;
}

std::string state_to_str(CELL_STATE st){
	switch(st){
	case ALIVE:
		return "ALIVE";
	case DEAD:
		return "DEAD";
	case DISA:
		return "DISA";
	case UNUSED:
		return "UNUSED";
	case FIX:
		return "FIX";
	case BLANK:
		return "BLANK";
	case DER:
		return "DER";
	case MUSUME:
		return "MUSUME";
	case AIR:
		return "AIR";
	case MEMB:
		return "MEMB";
	default:
		throw std::logic_error("Unknown state number:"+std::to_string((unsigned)st));
	}
}

int* __lat_x;
int* __lat_y;
int* __lat_z;
int* per_x_prev_idx;
int* per_x_next_idx;
int* per_y_prev_idx;
int* per_y_next_idx;
int* per_z_prev_idx;
int* per_z_next_idx;

void init_precalc_lat()
{
	static std::vector<int> latxv(pm->NX*3);
	static std::vector<int> latyv(pm->NY*3);
	static std::vector<int> latzv(pm->NZ*3);
    for (size_t i = 0; i < pm->NX * 3; i++) {
		latxv[i] = i%pm->NX;
	}
    for (size_t i = 0; i < pm->NY * 3; i++) {
		latyv[i] = i%pm->NY;
	}

    for (size_t i = 0; i < pm->NZ * 3; i++) {
		latzv[i] = i%pm->NZ;
	}
    __lat_x=&latxv[0];__lat_y=&latyv[0];__lat_z=&latzv[0];
}

void init_precalc_per()
{
	static std::vector<int> prex(pm->NX);static std::vector<int> nex(pm->NX);
	static std::vector<int> prey(pm->NY);static std::vector<int> ney(pm->NY);
	static std::vector<int> prez(pm->NZ);static std::vector<int> nez(pm->NZ);
    for (int k = 0; k<(int)pm->NX; k++) {
		int prev_x = k - 1;
		int next_x = k + 1;
		if (prev_x < 0) {
			prev_x += pm->NX;
		}
        if (next_x >= (int)pm->NX) {
			next_x -= pm->NX;
		}
		prex[k] = prev_x;
		nex[k] = next_x;
	}

    for (int k = 0; k<(int)pm->NY; k++) {
		int prev_y = k - 1;
		int next_y = k + 1;
		if (prev_y < 0) {
			prev_y += pm->NY;
		}
        if (next_y >= (int)pm->NY) {
			next_y -= pm->NY;
		}
		prey[k] = prev_y;
		ney[k] = next_y;
	}

    for (int l = 0; l < (int)pm->NZ; l++) {
		int prev_z = l-1, next_z = l+1;
		if (prev_z ==-1) {
			prev_z = 1;
		}
        if (next_z == (int)pm->NZ) {
			next_z = pm->NZ - 1;
		}
		prez[l] = prev_z;
		nez[l] = next_z;
	}

    per_x_prev_idx=&prex[0];per_x_next_idx=&nex[0];
    per_y_prev_idx=&prey[0];per_y_next_idx=&ney[0];
    per_z_prev_idx=&prez[0];per_z_next_idx=&nez[0];
}
