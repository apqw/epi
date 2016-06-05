#include "utils.h"
#include "define.h"
#include "cell.h"
#include <cmath>
double p_diff_x(double x1, double x2)
{
	using namespace cont;
	double diff = x1 - x2;
	if (fabs(diff) <= 0.5*LX)return diff;
	if (diff > 0)return diff - LX;
	if (diff < 0)return diff + LX;
	return diff;
}

double p_diff_y(double y1, double y2)
{
	using namespace cont;
	double diff = y1 - y2;
	if (fabs(diff) <= 0.5*LY)return diff;
	if (diff > 0)return diff - LY;
	if (diff < 0)return diff + LY;
	return diff;
}

double p_dist_sq(double x1, double y1, double z1, double x2, double y2, double z2)
{
	double diffx = p_diff_x(x1, x2);
	double diffy = p_diff_y(y1, y2);
	double diffz = z1 - z2;
	return diffx*diffx + diffy*diffy + diffz*diffz;
}

double p_cell_dist_sq(const Cell* c1, const Cell* c2)
{
	return p_dist_sq(c1->x(), c1->y(), c1->z(), c2->x(), c2->y(), c2->z());
}

double min0(double a) {
	return a > 0 ? a : 0;
}


bool no_double_count(Cell* c1, Cell* c2) {
	return c1->get_index() > c2->get_index();
}

int __lat_x[cont::NX*3];
int __lat_y[cont::NY * 3];
int __lat_z[cont::NZ * 3];
int per_x_prev_idx[cont::NX];
int per_x_next_idx[cont::NX];
int per_y_prev_idx[cont::NY];
int per_y_next_idx[cont::NY];

void init_precalc_lat()
{
	for (int i = 0; i < cont::NX * 3; i++) {
		__lat_x[i] = i%cont::NX;
	}
	for (int i = 0; i < cont::NY * 3; i++) {
		__lat_y[i] = i%cont::NY;
	}

	for (int i = 0; i < cont::NZ * 3; i++) {
		__lat_z[i] = i%cont::NZ;
	}
}

void init_precalc_per()
{
	using namespace cont;
	for (int k = 0; k<NX; k++) {
		int prev_x = k - 1;
		int next_x = k + 1;
		if (prev_x < 0) {
			prev_x += NX;
		}
		if (next_x >= NX) {
			next_x -= NX;
		}
		per_x_prev_idx[k] = prev_x;
		per_x_next_idx[k] = next_x;
	}

	for (int k = 0; k<NY; k++) {
		int prev_y = k - 1;
		int next_y = k + 1;
		if (prev_y < 0) {
			prev_y += NY;
		}
		if (next_y >= NY) {
			next_y -= NY;
		}
		per_y_prev_idx[k] = prev_y;
		per_y_next_idx[k] = next_y;
	}
}
