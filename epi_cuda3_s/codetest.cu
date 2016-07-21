#include "codetest.h"
#include "CellDataMan.h"
__global__ void dev_test(thrust::device_ptr<CellPos> cpos){

	cpos[0] = { 133343., 2., 3., 4. };


}

__global__ void periodic_diff_test(){
	printf("pdiff:%f - %f = %f old:%f\n", 20.0, threadIdx.x*0.112, p_diff_x(20.0, threadIdx.x*0.112), p_diff_x_old(20.0, threadIdx.x*0.112));
}
void codetest(){
	
	dbgexec([]{
		printf("dbg exec %d\n", 82179837928);

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
	printf("f %f\n", (real)pp(10, 10, 10));
	pp(10, 10, 10) = 0.0101010;
	printf("f %f\n", (real)pp(10, 10, 10));
	pp.fill(1919.0);
	printf("f %f\n", (real)pp(10, 10, 10));
	CArr3DMulti<real, 10, 10, 20, 3> mull;
	printf("f %f\n", (real)pp(10, 10, 10));
	CellConnectionData cchh;
	cchh.gj.init();
	cchh.gj.emplace(13234, 35435.55f);
	cchh.gj.emplace(13234, 35435.55f);
	cchh.gj.emplace(13234, 35435.55f);
	cchh.gj.emplace(13235, 35435.55f);
	cchh.gj.emplace(13236, 35435.55f);
	for (int i = 0; i < 1000; i++){
		cchh.gj.emplace(135 + i, 335435.55f);
	}
	printf("hashh %f\n", cchh.gj.at(703));
	cchh.gj.destroy();
	static_assert(std::is_trivially_copyable<IntegerHashmap<real>>::value, "real hashmap is not trivially copyable");
	/*
	for (int i = 0; i < 1000; i++){
	cchh.gj.remove(135 + i);
	}
	*/
	//pause();
	//((CellConnectionData)a.connection_data[0]).gj.emplace(135, 393894.4f);
	//a.check_initialized();
	periodic_diff_test << <1, 512 >> >();
	cudaDeviceSynchronize();
	//pause();
}

void env_check(){
	static_assert(std::is_standard_layout<IntegerHashmap<real>>::value, "IntegerHashmap<real> is not standard layout.\n");
	static_assert(std::is_trivially_copyable<IntegerHashmap<real>>::value, "IntegerHashmap<real> is not trivially copyable.\n");
}