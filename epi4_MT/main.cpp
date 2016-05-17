#include "CellData.h"
#include <stdio.h>
#include <cstdlib>
#include <cstring>
int main() {
	CellData* a = new CellData();
	std::memset(a->pos(), 0, sizeof(int)*MAX_CELL_NUM * 3);
	a->pos()[0][0] = 1919;
	//a->pos.old()[0][0] = 8181;
	a->pos.swap();
	printf("%lf\n", a->pos.old()[0][0]);
	system("pause");
}