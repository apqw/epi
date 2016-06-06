#include "utils.h"
#include "define.h"
#include "cell.h"
#include <cmath>



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
