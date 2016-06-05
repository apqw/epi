#include "Vec.h"

int main() {
	using namespace aml;
	Vec<double, 3> aaa = { 1,2,3 };
	Vec<double, 3> c = { 0.000001,0.000001,0.000001 };
	for (int i = 0; i < 10000000; i++) {
		aaa += c*aaa+aaa / aaa;
	}
	auto a1 = (c + aaa*aaa*aaa*aaa + aaa / c)[0];
	auto a2 = (aaa*aaa*aaa*aaa + aaa / c)[0];

	printf("%lf %lf %lf %lf %d \n", aaa[0], aaa[1],a1,a2,a1==a2);
	system("pause");
}