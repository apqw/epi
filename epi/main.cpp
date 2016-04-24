#include <iostream>
#include "define.h"

int main() {
	Vec4d a(1, 3, 5, 7);
	Vec4d b(5, 7, 11, 188);
	Vec4d c = m256dintmod(b, a);
	VSet4<Vec4d> pp;
	//VSet4d pp;
	EPI::Field_Data asss(1000, EPI::C::NX, EPI::C::NY, EPI::C::NZ);
	EPI::Field_Data ass2(1000, EPI::C::NX, EPI::C::NY, EPI::C::NZ);
	EPI::Ca2P::ca_dynamics(&asss, &ass2);
	pp.x = sqrt(a + b);
	
    std::cout << pp.x[0]<<","<< pp.x[1]<<","<< pp.x[2]<<","<< pp.x[3] << std::endl;
	system("pause");
    return 0;
}
