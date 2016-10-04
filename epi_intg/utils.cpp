#include "utils.h"
#include "global.h"
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
