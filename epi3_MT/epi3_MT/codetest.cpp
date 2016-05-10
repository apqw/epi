#include "codetest.h"
#include "define.h"
#include "component.h"
#include "Cell.h"
//#include "Field.h"
#include "primitive_func.h"
#include "general_func.h"
#include <cassert>

void Vec3Test()
{
	Vec3<int> vint = Vec3<int>({ 1, 2, 3 });
	Vec3<int> vint2 = Vec3<int>({ -2, -3, -4 });

	auto vintp = vint + vint2;
	auto s_vintm = 2 * vint;
	auto vint_sm = vint * 2;
	assert(vintp[0] == -1 && vintp[1] == -1 && vintp[2] == -1);
	assert(s_vintm[0] == 2 && s_vintm[1] == 4 && s_vintm[2] == 6);
	assert(vint_sm[0] == 2 && vint_sm[1] == 4 && vint_sm[2] == 6);
	int varr1[3] = { 2,4,6 }; std::array<int, 3> vstdarr1 = { 2,4,6 }; auto vvec1 = std::vector<int>({ 2,4,6 });
	assert(vint_sm == varr1); assert(vint_sm == vstdarr1); assert(vint_sm == vvec1);
	assert(vint_sm == Vec3<int>({ 2, 4, 6 }));
	assert(vint_sm == s_vintm);
	Vec3<double> vd({ 0.2,0.4,0.5 });
	auto vd1 = 0.5*vd;
	assert(vd1 == Vec3<double>({ 0.2*0.5,0.4*0.5,0.5*0.5 }));
	auto vd2 = vd / 0.5;
	assert(vd2 == Vec3<double>({ 0.2 / 0.5,0.4 / 0.5,0.5 / 0.5 }));
	printf("%lf\n", vd.normSq());
	assert(vd.normSq() == 0.2*0.2 + 0.4*0.4 + 0.5*0.5);
//	assert(Vec3<double>({ 0.5, 0.5, 0.5 }) == 0.5);
	vd1 = vd2;
	vd1 += vd2;
	printf("norm:%lf\n", sqrt(vd.normSq()));
	//vint_sm.norm();
	std::array<DV<double>, 3> aa = { 0.1,0.2,0.4 };
	auto p = Vec3<DV<double>>(aa);
	printf("Vec3Test end.\n");
}
class TestFunc {
public:
	void operator()(int& p) {
		p *= p;
	}
};
void arr_map_test()
{
	Arr3D<int> a = { {{1,2,3},{4,5,6},{7,8,9}} };
	arr_map<TestFunc,3>(a);
	assert(a[0][0][0] == 1 && a[0][0][1] == 4 && a[0][0][2] == 9);
	assert(a[0][1][0] == 16 && a[0][1][1] == 25 && a[0][1][2] == 36);
	assert(a[0][2][0] == 49 && a[0][2][1] == 64 && a[0][2][2] == 81);
	printf("arr_map_test end.\n");
}

void cell_test() {
	Cell a(ALIVE, { 1,2,3 }, 100);
	Cell b(ALIVE,  { -2,2,3 }, 100);
	assert(a.pos == Vec3<DV<double>>({ 1, 2, 3 }) && a.radius == 100);
//	a.radius = 0.0;
	assert(a.pos*b.pos == Vec3<DV<double>>({ -2, 4, 9 }) && a.radius == 100);
	a.pos+= a.pos*b.pos;
	assert(a.pos == Vec3<DV<double>>({ 1, 2, 3 }));
	a.pos.foreach([](DV<double>& b) {b.update(); });
	assert(a.pos == Vec3<DV<double>>({ -1, 6, 12 }));
	printf("cell_test end.\n");
}

void conv_state_num_test() {
	/*
	assert(conv_state_num(ALIVE) == ALIVE_N);
	assert(conv_state_num(DEAD) == DEAD_N);
	assert(conv_state_num(DISA) == DISA_N);
	assert(conv_state_num(UNUSED) == UNUSED_N);
	assert(conv_state_num(FIX) == FIX_N);
	assert(conv_state_num(BLANK) == BLANK_N);
	assert(conv_state_num(DER) == DER_N);
	assert(conv_state_num(MUSUME) == MUSUME_N);
	assert(conv_state_num(AIR) == AIR_N);
	assert(conv_state_num(MEMB) == MEMB_N);
	*/
	printf("conv_state_num_test end.\n");
}

void cell_man_test() {
//	CellMan<CellRef> a(100);
//	a.add(CellRef(new Cell(ALIVE)));
	CellMan a; LWCellMan<400> lw;
	a.add_queue(std::make_shared<Cell>(ALIVE, std::initializer_list<DV<double>> {-1.1,-4.1,-9.1} , 100));
	a.add_queue(std::make_shared<Cell>(DISA));
	a.update();
	a.foreach([](CellPtr& c, int i) {
		c->pos += {1, 2, 3};
	
	});


	a.foreach([](CellPtr& c, int i) {
		printf("%lf %lf %lf\n", (double)c->pos[0], (double)c->pos[1], (double)c->pos[2]);
		c->update();
	});

	a.foreach([](CellPtr& c, int i) {
		printf("%lf %lf %lf\n", (double)c->pos[0], (double)c->pos[1], (double)c->pos[2]);
		//c->update();
	});

	a.foreach([&lw](CellPtr& c, int i) {
		lw.add(c);
	});

	lw.foreach([](Cell* c) {
		printf("%lf %lf %lf\n", (double)c->pos[0], (double)c->pos[1], (double)c->pos[2]);
		c->pos *= {1, 2, 3};
		c->update();
	});

	lw.foreach([](Cell* c) {
		printf("%lf %lf %lf\n", (double)c->pos[0], (double)c->pos[1], (double)c->pos[2]);
		c->pos *= {1, 2, 3};
		c->update();
	});

	a.remove_if([](CellPtr& c) {return c->pos[0] < 0.0; });
	//a.update();
	a.foreach([](CellPtr& c, int i) {
		printf("%lf %lf %lf\n", (double)c->pos[0], (double)c->pos[1], (double)c->pos[2]);
		//c->update();
	});

}

void list_man_test() {
	//ListMan<CellRef> a(100);
//
}

void CAS_test() {
	CAS<double> a = 0.0;
	double a2 = 0.0;
    std::atomic<double> a3 (0.0);
	DV<double> ap = 0.0;
	constexpr int loop = 81818812;
	printf("CAS testing...\n");
//#pragma omp parallel for
	for (int i = 0; i < loop; i++) {
		a += 1.0;
		a2 += 1.0;
		a3.store(a3.load() + 1);
		ap += 1.0;
	}
	assert(a == (double)loop);
	ap.update();
    printf("omp normal add:%lf CAS add:%lf atm add:%lf dv add:%lf\n", a2,a.v.load(std::memory_order_relaxed),a3.load(),(double)ap);
}

void init_test(std::ifstream & strm)
{

}
