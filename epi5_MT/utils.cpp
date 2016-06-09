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
int per_z_prev_idx[cont::NZ];
int per_z_next_idx[cont::NZ];

void init_precalc_lat()
{
    for (size_t i = 0; i < cont::NX * 3; i++) {
		__lat_x[i] = i%cont::NX;
	}
    for (size_t i = 0; i < cont::NY * 3; i++) {
		__lat_y[i] = i%cont::NY;
	}

    for (size_t i = 0; i < cont::NZ * 3; i++) {
		__lat_z[i] = i%cont::NZ;
	}
}

void init_precalc_per()
{
	using namespace cont;
    for (int k = 0; k<(int)NX; k++) {
		int prev_x = k - 1;
		int next_x = k + 1;
		if (prev_x < 0) {
			prev_x += NX;
		}
        if (next_x >= (int)NX) {
			next_x -= NX;
		}
		per_x_prev_idx[k] = prev_x;
		per_x_next_idx[k] = next_x;
	}

    for (int k = 0; k<(int)NY; k++) {
		int prev_y = k - 1;
		int next_y = k + 1;
		if (prev_y < 0) {
			prev_y += NY;
		}
        if (next_y >= (int)NY) {
			next_y -= NY;
		}
		per_y_prev_idx[k] = prev_y;
		per_y_next_idx[k] = next_y;
	}

    for (int l = 0; l < (int)NZ; l++) {
		int prev_z = l-1, next_z = l+1;
		if (prev_z ==-1) {
			prev_z = 1;
		}
        if (next_z == (int)NZ) {
			next_z = NZ - 1;
		}
		per_z_prev_idx[l] = prev_z;
		per_z_next_idx[l] = next_z;
	}
}
