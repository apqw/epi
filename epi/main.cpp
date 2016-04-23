#include <iostream>
#include "define.h"

int main() {
	Vec4d a(1, 3, 5, 7);
	Vec4d b(5, 7, 11, 188);
	Vec4d c = m256dintmod(b, a);
	VSet4<Vec4d> pp;
	//VSet4d pp;
	pp.x = sqrt(a + b);
	
    std::cout << pp.x[0]<<","<< pp.x[1]<<","<< pp.x[2]<<","<< pp.x[3] << std::endl;
	system("pause");
    return 0;
}
