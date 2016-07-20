#include "codetest.h"
#include "CellDataMan.h"
__global__ void dev_test(thrust::device_ptr<CellPos> cpos){

	cpos[0] = { 133343., 2., 3., 4. };


}
void codetest(){
	dbgexec([]{
		printf("dbg exec %d\n",82179837928);

	});
	CellDataMan(a);
	printf("%f\n", ((CellPos)a.pos.current()[0]).x);
	dev_test << <1, 1 >> >(a.pos.current());
	printf("%f\n", ((CellPos)a.pos.current()[0]).x);
	a.cell_phase_switch();
	printf("%f\n", ((CellPos)a.pos.current()[0]).x);
	a.cell_phase_switch();
	printf("%f\n", ((CellPos)a.pos.current()[0]).x);
	Field3D<real> pp;
	printf("f %f\n", (real)pp(10,10,10));
	pp(10, 10, 10) = 0.0101010;
	printf("f %f\n", (real)pp(10, 10, 10));
	pp.fill(1919.0);
	printf("f %f\n", (real)pp(10, 10, 10));
	CArr3DMulti<real, 10, 10, 20, 3> mull;
	printf("f %f\n", (real)pp(10, 10, 10));
	//a.check_initialized();
}