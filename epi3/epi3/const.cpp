#include "const.h"
constexpr double constabs(double a) {
	return a > 0 ? a : -a;
}
void static_error_check() {
	static_assert(constabs(dx - dh) <= EPS && constabs(dy - dh) <= EPS && constabs(dz - dh) <= EPS,
		"error:grid size must be near dh * dh * dh\n");

	constexpr double tdiv = agki_max * (1.0 - stoch_div_time_ratio) / (eps_kb * u0);
	constexpr double tseparate = 2 * R_max / eps_L;
	static_assert(tdiv >= tseparate, "error: division period must be greater than the time for the completion of cell division\n");
}